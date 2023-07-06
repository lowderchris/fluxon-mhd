
from sunpy.net import attr, attrs, Fido
from sunpy.net.dataretriever import GenericClient
from sunpy.net.dataretriever.client import QueryResponse
from sunpy.net.scraper import Scraper
from sunpy.time import TimeRange
import os

class PUNCHType(attr.SimpleAttr):
    """
    MPM: Magnetic Pseudo-Moments
    """
    pass

class PUNCHLevel(attr.SimpleAttr):
    """
    Preliminary synthetic PUNCH data generated from GAMERA model data.

    level3 - Level 3 total brightness and polarized brightness at 24-minute cadence
    level3_trefoil - Quasi level 3 total brightness and polarized brightness, with the trefoil FOV at 8-minute cadence
    quickPUNCH - Lower resolution quick products with total brightness only

    Note that the level3_trefoil and quickPUNCH products are at an 8-minute cadence for purposes of illustration, and contain a CME of unusual speed. 
    The level3 products are synthetic assmblies of three sets of spacecraft rotation trefoil patterns, and are a more accurate representation.
    """
    pass

class PUNCHClient(GenericClient):

    top_lvl = 'https://data.boulder.swri.edu/lowder/PUNCH/synthetic_data/{PUNCHLevel}/'
    baseurl = top_lvl + "PUNCH_L(\d){1}_(\w){3}_(\d){14}.fits"
    pattern = top_lvl + "PUNCH_L{Level:1d}_{PUNCHType:3l}_{year:4d}{month:2d}{day:2d}{hour:2d}{minute:2d}{second:2d}.fits"
    info_url = top_lvl


    @classmethod
    def register_values(cls):
        from sunpy.net import attrs as a

        adict = {attrs.Instrument: [('PUNCH', 'Polarimeter to UNify the Corona and Heliosphere')],
                attrs.Source: [('SwRI', 'Southwest Research Institute')],
                attrs.Provider: [('SwRI', 'Southwest Research Institute')],
                attrs.Level: [('2', "level2"), ('3', 'level3')],
                PUNCHType: [('MPM', 'Magnetic Pseudo-Moments'), ('NQM', 'Non-Quadrupolar Moments'), ('WQM', 'Wavelet Moments')], #This is almost certainly complete nonsense
                PUNCHLevel: [('level3', "Regular Level 3"),('level3_trefoil','The Trefoil Image'), ("quickPUNCH", "A smaller quicklook image")]
                }
        return adict
    
    @classmethod
    def _can_handle_query(cls, *query):
        required = {attrs.Instrument, attrs.Time}

        optional = {attrs.Level, attrs.Source, attrs.Provider, PUNCHType, PUNCHLevel}

        all_attrs = {type(x) for x in query}
        return required.issubset(all_attrs) and all_attrs.issubset(required.union(optional))

    def search(self, *args, **kwargs):
        """
        Query this client for a list of results.

        Parameters
        ----------
        \\*args: `tuple`
            `sunpy.net.attrs` objects representing the query.
        \\*\\*kwargs: `dict`
             Any extra keywords to refine the search.

        Returns
        -------
        A `QueryResponse` instance containing the query result.
        """

        baseurl, pattern, matchdict = self.pre_search_hook(*args, **kwargs)
        tr = TimeRange(matchdict['Start Time'], matchdict['End Time'])
        filesmeta_list = []
        for key in matchdict:
            try: # if the value is a list, we need to iterate over it
                for value in matchdict[key]:
                    kwargs_dict = {key: value}
                    this_baseurl = baseurl.replace("{{{}}}".format(key), value)
                    scraper = Scraper(this_baseurl, regex=True, **kwargs_dict)
                    filesmeta_list.extend(scraper._extract_files_meta(tr, extractor=pattern, matcher=matchdict))
            except TypeError:
                pass

        filesmeta = sorted(filesmeta_list, key=lambda k: k['url'])
        return QueryResponse(filesmeta, client=self)


def create_punchdir():
    """
    Create a directory to store PUNCH data
    """
    save_dir = os.getcwd().replace('mhd', 'data') + "/PUNCH_DL"
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

#run code if this is the main file
if __name__ == '__main__':

    create_punchdir()

    date_start = '2023-07-04T00:00:00'
    date_end = '2023-7-04T01:00:00'

    res = Fido.search(attrs.Instrument('punch'), attrs.Time(date_start, date_end), attrs.Level('3')) #, PUNCHLevel('level3_trefoil'))
    print(res)
    # print(save_dir)
    # out = Fido.fetch(res, path=save_dir)
