import pandas as pd
import numpy as np

def merge_reduced(result_df: pd.DataFrame,
                  df: pd.DataFrame,
                  map_cols: dict[str, str],
                  reduction: str,
                  merge_on: list[str] | str=['ImageNumber', 'ObjectNumber'],
                  how: str='inner',
                  ):
    '''Merge df into result_df, reducing as needed.
    Only columns in map_cols are considered, and are renamed'''
    return result_df.merge(
        df.groupby(merge_on)[list(map_cols.keys())].aggregate(reduction).reset_index().rename(columns=map_cols),
        on=merge_on,
        how=how,
    )

def merge_result(result_df: pd.DataFrame,
                 df: pd.DataFrame,
                  map_cols: dict[str, str],
                  merge_on: list[str] | str=['ImageNumber', 'ObjectNumber'],
                  how: str='inner',
                  ):
    return result_df.merge(
        df[merge_on + list(map_cols.keys())].rename(columns=map_cols),
        on=merge_on,
        how=how,
    )


def analyze(
    datafile,
    parsers=[],
    extra_columns=[],
    previous_result=None,
    merge_fcn=None,
    region='',
    reduce=False,
    ):
    for parser in parsers:
        parser.peek_file(datafile)
    data = pd.read_csv(datafile,
                       usecols=sum(
                        (parser.get_columns() for parser in parsers),
                        extra_columns,
                       )
                       )

    extra_data = []
    for parser in parsers:
        extra = parser.analyze(data)
        if extra is not None:
            extra_data.append(extra)

    if previous_result is not None and merge_fcn:
        for parser in parsers:
            if reduce:
                previous_result = parser.merge_reduced_result(
                    previous_result, data, region, merge_fcn)

            else:
                previous_result = parser.merge_result(
                    previous_result, data, region, merge_fcn)
        return previous_result, extra_data

    return data, extra_data



class CellProfilerParser():
    def peek_file(self, datafile):
        # optionally open and inspect file if e.g. needed to get columns
        pass

    def get_columns(self) -> list[str]:
        return []

    def analyze(self, df) -> pd.DataFrame | None:
        '''Add derived columns to df in place, return new df if needed.'''
        pass

    def merge_result(self, result, df, region, merge_fcn):
        '''Merge df with result for the given region.
        merge_fcn will be called with result, df, map_cols.'''
        pass

    def merge_reduced_result(self, result, df, region, merge_fcn):
        '''Merge df with result for the given region.
        merge_fcn will be called with result, df, map_cols, reduction.'''
        pass


class RDFParser(CellProfilerParser):
    def __init__(self, id_vars=[]):
        self.columns = []
        self.id_vars = id_vars

    def peek_file(self, datafile):
        self.columns = [
            column for column in open(datafile).readline().strip().split(',')
                if column.startswith('RDF_')
        ]

    def get_columns(self) -> list[str]:
        return self.columns + self.id_vars

    def analyze(self, df) -> pd.DataFrame | None:
        rdf_raw = df.melt(id_vars=self.id_vars)

        intensity = rdf_raw[rdf_raw.variable.str.startswith('RDF_Intensity')].reset_index(drop=True)
        extract = intensity.variable.str.extract(r'RDF_Intensity_C(\d)_R([-0-9]+)')
        intensity = intensity.assign(
            channel=extract[0].astype(int),
            radius=extract[1].astype(int),
        ).rename(columns={'value': 'intensity'}).drop(columns='variable')

        counts = rdf_raw[rdf_raw.variable.str.startswith('RDF_Count')].reset_index(drop=True)
        extract = counts.variable.str.extract(r'RDF_Counts_R([-0-9]+)')

        counts = counts.assign(
            radius=extract[0].astype(int)
        ).rename(columns={'value': 'counts'}).drop(columns='variable')

        return intensity.merge(counts, on=self.id_vars + ['radius'])


    def merge_result(self, result, df, region, merge_fcn):
        return result

    def merge_reduced_result(self, result, df, region, merge_fcn):
        return result


class CorrelationParser(CellProfilerParser):
    def __init__(self, measures, reduce=False):
        self.columns = []
        self.measures = measures
        self.reduce = reduce

    def peek_file(self, datafile):
        measures = tuple(f'Correlation_{measure}' for measure in self.measures)
        self.columns = [
            column for column in open(datafile).readline().strip().split(',')
            if column.startswith(measures)
        ]

    def get_columns(self) -> list[str]:
        return self.columns + (['AreaShape_Area'] if self.reduce else [])

    def merge_result(self, result, df, region, merge_fcn):
        map_cols = {
            column: f'{region}_{column[12:]}'
            for column in self.columns
        }
        return merge_fcn(result, df, map_cols)

    def merge_reduced_result(self, result, df, region, merge_fcn):
        map_cols = {
            column: f'Mean_{region}_{column[12:]}'
            for column in self.columns
        }
        map_cols['AreaShape_Area'] = '_tmp_total_area'
        result = merge_fcn(result, df, map_cols, 'sum')
        result[list(map_cols.values())] /= result['_tmp_total_area'].to_numpy()[:, None]
        result = result.drop(columns='_tmp_total_area')

        return result


class RimEnrichmentParser(CellProfilerParser):
    def __init__(self, images: list, area_normalization, bins, total_bins):
        self.images = images
        self.area_normalization = area_normalization
        self.bins = bins
        self.total_bins = total_bins


    def _bins(self):
        return range(self.total_bins, self.total_bins-self.bins, -1)

    def get_columns(self) -> list[str]:
        return [
            f'RadialDistribution_FracAtD_{image}_{bin}of{self.total_bins}'
            for image in (self.images + [self.area_normalization])
            for bin in self._bins()
        ]

    def analyze(self, df) -> pd.DataFrame | None:
        relative_area = df[[
            f'RadialDistribution_FracAtD_{self.area_normalization}'
            f'_{bin}of{self.total_bins}'
                for bin in self._bins()
        ]].sum(axis=1)

        for image in self.images:
            df[f'{image}_Rim_Enrichment'] = df[[
                f'RadialDistribution_FracAtD_{image}_{bin}of{self.total_bins}'
                for bin in self._bins()
            ]].sum(axis=1) / relative_area

    def merge_result(self, result, df, region, merge_fcn):
        '''Merge df with result for the given region.
        merge_fcn will be called with result, df, map_cols.'''
        map_cols = {
            f'{image}_Rim_Enrichment': f'{region}_{image}_Rim_Enrichment'
            for image in self.images
        }

        return merge_fcn(result, df, map_cols)


    def merge_reduced_result(self, result, df, region, merge_fcn):
        '''Merge df with result for the given region.
        merge_fcn will be called with result, df, map_cols, reduction.'''
        map_cols = {
            f'{image}_Rim_Enrichment': f'Mean_{region}_{image}_Rim_Enrichment'
            for image in self.images
        }

        return merge_fcn(result, df, map_cols, 'mean')


class IntensityParser(CellProfilerParser):
    def __init__(self, measures=['Mean'], images=list(), locations=list(), total_intens=False):
        self.measures = measures
        self.images = images
        self.locations = locations
        self.total_intens = total_intens

    def get_columns(self) -> list[str]:
        result = []

        if self.total_intens:
            result = ['AreaShape_Area']
            self.measures.append('Mean')

        result += [
            f'Intensity_{measure}Intensity_{image}'
                for measure in self.measures
                for image in self.images
        ]
        if self.locations:
            result += [
                f'Location_CenterMassIntensity_{coord}_{image}'
                for coord in ('X', 'Y')
                for image in self.images
            ]

        return result

    def analyze(self, df) -> pd.DataFrame | None:
        '''Add derived columns to df in place, return new df if needed.'''
        if self.total_intens:
            for image in self.images:
                df[f'Intensity_TotalIntensity_{image}'] = (
                    df[f'Intensity_MeanIntensity_{image}']
                        * df['AreaShape_Area'])

    def _map_cols(self, df, prefix):
        map_cols = {
            column: column.replace('Intensity', prefix, 1)
            for column in df.columns
            if column.startswith('Intensity')
        }

        if self.locations:
            map_cols.update({
                column: column.replace('Location', prefix)
                for column in df.columns
                if column.startswith('Location')
            })

        return map_cols


    def merge_result(self, result, df, region, merge_fcn):
        '''Merge df with result for the given region.
        merge_fcn will be called with result, df, map_cols.'''
        map_cols = self._map_cols(df, region)
        return merge_fcn(result, df, map_cols)

    def merge_reduced_result(self, result, df, region, merge_fcn):
        '''Merge df with result for the given region.
        merge_fcn will be called with result, df, map_cols, reduction.'''

        map_cols = self._map_cols(df, f'Mean_{region}')

        result = merge_fcn(result, df, map_cols, 'mean')

        if self.total_intens:
            map_cols = {
                'AreaShape_Area': '_tmp_total_area',
            }
            map_cols.update({
                f'Intensity_TotalIntensity_{image}': f'Total_{region}_TotalIntensity_{image}'
                for image in self.images
            })

            result = merge_fcn(result, df, map_cols, 'sum')

            for image in self.images:
                result[f'MeanWeighted_{region}_Intensity_{image}'] = (
                    result[f'Total_{region}_TotalIntensity_{image}'] / 
                    result['_tmp_total_area']
                )

            result = result.drop(columns='_tmp_total_area')

        return result

        

class BlankParser(CellProfilerParser):
    def __init__(self, columns=list()):
        self.columns = columns

    def get_columns(self):
        return self.columns

    '''Pass through any columns, useful for joining datasets'''
    def merge_result(self, result, df, region, merge_fcn):
        '''Merges all columns in df.'''
        return merge_fcn(result, df, dict(zip(self.columns, self.columns)))

    def merge_reduced_result(self, result, df, region, merge_fcn):
        raise ValueError('Operation not allowed for BlankParser')


class CountingParser(CellProfilerParser):
    '''Count number of observations for merging.'''
    def get_columns(self) -> list[str]:
        return ['ObjectNumber']

    def merge_result(self, result, df, region, merge_fcn):
        '''Do nothing on non-reduced merge.'''
        raise ValueError('Operation not allowed for CountingParser')

    def merge_reduced_result(self, result, df, region, merge_fcn):
        map_cols = {'ObjectNumber': f'Count_{region}'}
        return merge_fcn(result, df, map_cols, 'count')


class ImageParser(CellProfilerParser):
    def __init__(self, regex=None, debug_regex=False):
        self.regex = regex
        self.map_cols = []
        self.debug = debug_regex

    def get_columns(self) -> list[str]:
        return ['Metadata_FileLocation', 'Metadata_Series', 'ImageNumber']

    def analyze(self, df) -> pd.DataFrame | None:
        if self.regex is None:
            return
        parsed = df['Metadata_FileLocation'].str.extract(self.regex)

        if self.debug and parsed.isna().any(axis=None):
            print(df.loc[parsed.isna().any(axis=1), 'Metadata_FileLocation'].unique())
            raise ValueError("Regex failed")

        self.map_cols = {column: column for column in parsed.columns}
        self.map_cols['Metadata_FileLocation'] = 'Metadata_FileLocation'

        # joining doesn't modify the df in place, use columns instead
        for col in parsed.columns:
            df[col] = parsed[col]

    def merge_result(self, result, df, region, merge_fcn):
        return merge_fcn(result, df, self.map_cols)

    def merge_reduced_result(self, result, df, region, merge_fcn):
        raise ValueError('Operation not allowed for FilelocationParser')


class ShapeParser(CellProfilerParser):
    def get_columns(self):
        return ['AreaShape_Area',
                'AreaShape_Eccentricity',
                'AreaShape_Perimeter',
                ]

    def analyze(self, df):
        df['AreaShape_Circularity'] = (
            4 * np.pi * df['AreaShape_Area'] /
                df['AreaShape_Perimeter'] ** 2
        )

    def merge_result(self, result, df, region, merge_fcn):
        columns = self.get_columns() + ['AreaShape_Circularity']
        map_cols = {
            column: column.replace('AreaShape', region)
            for column in columns
        }

        return merge_fcn(result, df, map_cols)

    def merge_reduced_result(self, result, df, region, merge_fcn):
        map_cols = {
            f'AreaShape_{measure}': f'Mean_{region}_{measure}'
            for measure in ('Area', 'Circularity', 'Perimeter', 'Eccentricity')
        }
        result = merge_fcn(result, df, map_cols, 'mean')

        map_cols = {
            'AreaShape_Area': f'Total_{region}_Area',
        }
        result = merge_fcn(result, df, map_cols, 'sum')

        return result


# Intensity gives center of mass
