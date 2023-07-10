
from sunpy.net import attr, attrs, Fido
a = attrs
from sunpy.net.dataretriever import GenericClient
from sunpy.net.dataretriever.client import QueryResponse
from sunpy.net.scraper import Scraper
from sunpy.time import TimeRange
import os
from sunpy.net.attr import SimpleAttr

## Custom attributes for this dataset ##
class PUNCHSeries(attr.SimpleAttr):
    """
    NQ is NFI quickpunch, 
    WQ is the combined WFI+NFI quickpunch, 
    MP is the polarized mosaic product. 

    (NQ, WQ, MP)

    """

class PUNCHClient(GenericClient):

    top_lvl = 'https://data.boulder.swri.edu/lowder/PUNCH/synthetic_data/'
    baseurl = top_lvl + "{ExtentType}/PUNCH_L(\d){1}_(\w){3}_(\d){14}.fits"
    pattern = top_lvl + "{ExtentType}/PUNCH_L{Level:1d}_{PUNCHSeries:2l}{Detector:1l}_{year:4d}{month:2d}{day:2d}{hour:2d}{minute:2d}{second:2d}.fits"
    info_url = top_lvl

    @classmethod
    def register_values(cls):

        # For every attribute, assign a list of all possible (value, description) tuples
        # Should I rename any of the custom attributes?
        adict = {
                attrs.Source: [('PUNCH', 'Polarimeter to UNify the Corona and Heliosphere')],
                attrs.Provider: [('SwRI', 'Southwest Research Institute')], # should this be a different thing than "Source"?


                attrs.Level: [('0', "level0: Raw from camera"),
                              ('1', "level1: *Description*"), 
                              ('2', "level2: *Description*"), 
                              ('3', 'level3: *Description*"')], 
                
                attrs.ExtentType: [('level3', "Regular Level 3"), 
                                   ('level3_trefoil','The Trefoil Image'), 
                                   ("quickPUNCH", "A smaller quicklook image")], 

                #     Preliminary synthetic PUNCH data generated from GAMERA model data.

                #     level3 - Level 3 total brightness and polarized brightness at 24-minute cadence
                #     level3_trefoil - Quasi level 3 total brightness and polarized brightness, with the trefoil FOV at 8-minute cadence
                #     quickPUNCH - Lower resolution quick products with total brightness only

                #     Note that the level3_trefoil and quickPUNCH products are at an 8-minute cadence for purposes of illustration, and contain a CME of unusual speed. 
                #     The level3 products are synthetic assmblies of three sets of spacecraft rotation trefoil patterns, and are a more accurate representation.
                #     (level3, level3_trefoil, quickPunch)

                # NQ is NFI quickpunch, WQ is the combined WFI+NFI quickpunch, MP is the polarized mosaic product. 
                PUNCHSeries: [('MP', 'Polarized Mosaic'), 
                              ('NQ', 'NFI quickpunch'), 
                              ('WQ', 'WFI + NFI quickpunch')], 

                ## The third character is the observatory designation, a number 1-4 for the four spacecraft, or N for a NFI-only mosaic, or M for a mosaic product
                attrs.Detector: [('M', "PUNCH Mosaic"), ('N', "NFI from PUNCH"), ('W', "WFI from PUNCH"), 
                                 ('1', "WFI 1"), ('2', "WFI 2"), ('3', "WFI 3")], 

                # attrs.Instrument: [('WFI', 'Wide Field Imager'), ('NFI', 'Narrow Field Imager'), ('M', 'Mosaic')],

                }
        return adict
    
    @classmethod
    def _can_handle_query(cls, *query):
        required = {attrs.Time}

        optional = {attrs.Source, attrs.Provider, attrs.Level, attrs.ExtentType, PUNCHSeries, attrs.Detector} #, attrs.Instrument}

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

        ## Original, with verb=True
        # verb=False
        # baseurl, pattern, matchdict = self.pre_search_hook(*args, **kwargs)
        # if verb: print("Args: ", args)
        # if verb: print("Baseurl: ", baseurl)
        # # for ii, key in enumerate(matchdict.keys()):
        # #     print(key," - ", matchdict[key])
        # #     print(args[ii]," - ", args[ii].value)
        # if verb: print("Matchdict: ", matchdict)
        # # if verb: print("Matchdict: ", matchdict["Instrument"][0])
        # # if verb: print("SearchArg: ", args[0].value)
        # scraper = Scraper(baseurl, regex=True)
        # # import pdb; pdb.set_trace()
        # tr = TimeRange(matchdict['Start Time'], matchdict['End Time'])
        # filesmeta = scraper._extract_files_meta(tr, extractor=pattern,
        #                                         matcher=matchdict)
        # filesmeta = sorted(filesmeta, key=lambda k: k['url'])
        # metalist = []
        # for i in filesmeta:
        #     rowdict = self.post_search_hook(i, matchdict)
        #     metalist.append(rowdict)
        # # import pdb; pdb.set_trace()
        # return QueryResponse(metalist, client=self)
    
        ## Custom. This lets the url have variables in it, defined as {variable_name} in the url
        baseurl, pattern, matchdict = self.pre_search_hook(*args, **kwargs)
        tr = TimeRange(matchdict['Start Time'], matchdict['End Time'])
        filesmeta_list = []
        for key in matchdict:
            try: # if the value is a list, we need to iterate over it
                for value in matchdict[key]:
                    if key in baseurl:
                        this_baseurl = baseurl.replace("{{{}}}".format(key), value) ############
                        scraper = Scraper(this_baseurl, regex=True) ############
                        filesmeta_list.extend(scraper._extract_files_meta(tr, extractor=pattern, matcher=matchdict)) ############
            except TypeError:
                pass

        filesmeta = sorted(filesmeta_list, key=lambda k: k['url'])
        return QueryResponse(filesmeta, client=self)
    

    @classmethod
    def _get_match_dict(cls, *args, **kwargs):
        """
        Constructs a dictionary using the query and registered Attrs that represents
        all possible values of the extracted metadata for files that matches the query.
        The returned dictionary is used to validate the metadata of searched files
        in :func:`~sunpy.net.scraper.Scraper._extract_files_meta`.

        Parameters
        ----------
        \\*args: `tuple`
            `sunpy.net.attrs` objects representing the query.
        \\*\\*kwargs: `dict`
             Any extra keywords to refine the search.

        Returns
        -------
        matchdict: `dict`
            A dictionary having a `list` of all possible Attr values
            corresponding to an Attr.
        """
        regattrs_dict = cls.register_values()
        matchdict = {}
        for i in regattrs_dict.keys():
            attrname = i.__name__
            # only Attr values that are subclas of Simple Attr are stored as list in matchdict
            # since complex attrs like Range can't be compared with string matching.
            if issubclass(i, SimpleAttr):
                matchdict[attrname] = []
                for val, _ in regattrs_dict[i]:
                    # print(val)
                    matchdict[attrname].append(val)
        # print("\n")
        for elem in args:
            if isinstance(elem, a.Time):
                matchdict['Start Time'] = elem.start
                matchdict['End Time'] = elem.end
            elif hasattr(elem, 'value'):
                # This is the only custom line, which lets the search be case-insensitive
                matchdict[elem.__class__.__name__] = [str(elem.value).lower(), str(elem.value).upper(), str(elem.value)] #############
            elif isinstance(elem, a.Wavelength):
                matchdict['Wavelength'] = elem
            else:
                raise ValueError(
                    "GenericClient can not add {} to the rowdict dictionary to"
                    "pass to the Client.".format(elem.__class__.__name__))
        return matchdict


## ## Helper Functions #################################################

def create_punchdir(basedir=None):
    """
    Create a directory to store PUNCH data
    """
    if basedir is None:
        basedir = os.getcwd().replace('mhd', 'data') 
    save_dir = basedir + "/PUNCH_COMP/"
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    return save_dir + "{extenttype}" # has to be all lowercase


def get_all_files(directory):
    import glob
    pattern = directory + '/**/*'  # Recursive pattern to match all files
    file_list = glob.glob(pattern, recursive=True)
    return file_list


def unpack_punchfile(filename, save_dir=None):
    """
    Unpacks a PUNCH file and saves it to a directory
    """
    if save_dir is None:
        # save_dir = filename.replace('.fits', '_open.fits')
        save_dir = filename.replace('PUNCH_COMP', 'PUNCH_OPEN')
        if not os.path.exists(os.path.dirname(save_dir)):
            os.makedirs(os.path.dirname(save_dir))


    import astropy.io.fits as fits
    import sunpy.map
    
    with fits.open(filename) as hdul:
        this_map = sunpy.map.Map(hdul[1].data, hdul[1].header)
        this_map.save(save_dir, filetype='fits', overwrite=True)


def unpack_all_punchfiles(directory, delete=True):
    #get the filenames of the files that were going to be downloaded but are already in the directory
    print("\n\nUnpacking PUNCH files...")
    paths = get_all_files(os.path.dirname(directory))
    fits_paths = [x for x in paths if x.endswith('.fits')]
    needed_paths = [x for x in fits_paths if not '_open' in x]
    for file in needed_paths:
        unpack_punchfile(file)
        if delete:
            os.remove(file)
    if delete:
        os.rmdir(os.path.dirname(file))
        os.rmdir(os.path.dirname(directory))
    print("Done!")

#######################################################################################################################################


def test_PUNCHClient(date_start = '2023-07-04T00:00:00', date_end = '2023-7-04T01:00:00', unpack=True, delete=False):
    """Downloads an example PUNCH file and saves it to a directory

    Args:
        date_start (str, optional): _description_. Defaults to '2023-07-04T00:00:00'.
        date_end (str, optional): _description_. Defaults to '2023-7-04T01:00:00'.
    """
    save_dir = create_punchdir()

    #{attrs.Time, attrs.Source, attrs.Provider, attrs.Level, attrs.ExtentType, PUNCHSeries, attrs.Detector}

    res = Fido.search(attrs.Source('punch'), attrs.Time(date_start, date_end), attrs.ExtentType('level3'))

    print(res)    
    user_input = input(f"Save directory: {save_dir}\n\nDo you want to proceed to download? [Y]/N/E[xit]: ").strip().lower()
    if user_input == 'y' or user_input == '':
        out = Fido.fetch(res, path=save_dir, overwrite=False)
    elif user_input in ['e', 'exit', 'q', 'quit']:
        print("Exiting program.\n\n")
        return None
    else:
        out = None
        print("\t\t\t\t\tDownload cancelled.", end="")

    if unpack:
        unpack_all_punchfiles(save_dir, delete=delete)    

    print("\n\nProgram complete.\n\n\n\n")
    return out
    


#run code if this is the main file
if __name__ == '__main__':
    test_PUNCHClient(unpack=True)



        

