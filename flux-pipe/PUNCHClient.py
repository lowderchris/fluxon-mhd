
import sunpy.net
from sunpy.net import attr, attrs
from sunpy.net.dataretriever import GenericClient


class PUNCHType(attr.SimpleAttr):
    """
    MPM: Magnetic Pseudo-Moments
    """
    pass

class PUNCHLevel(attr.SimpleAttr):
    """
    
    """
    pass


class PUNCHClient(GenericClient):


    folder = "/{PUNCHLevel}/" #"/level3/"  #"level3"  # "level3" #
    # print(PUNCHLevel)
    #Model url:   https://data.boulder.swri.edu/lowder/PUNCH/synthetic_data/level3/PUNCH_L3_MPM_20230704151200.fits
    #Model url:   https://data.boulder.swri.edu/lowder/PUNCH/synthetic_data/level3_trefoil/PUNCH_L3_MPM_20230704151200.fits
    top_lvl = 'https://data.boulder.swri.edu/lowder/PUNCH/synthetic_data'
    baseurl = top_lvl + folder + "PUNCH_L3_MPM_(\d){14}.fits"
    # print(baseurl)
    # import pdb; pdb.set_trace()
    # baseurl = baseurl.replace("{PUNCHLevel}", "level3")
    # print(baseurl)
    pattern = top_lvl + "/{PUNCHLevel}/PUNCH_L{Level:1d}_{PUNCHType:3l}_{year:4d}{month:2d}{day:2d}{hour:2d}{minute:2d}{second:2d}.fits"

    info_url = baseurl
    __name__ = 'PUNCHClient'

    # # @classmethod
    # def name_it(cls, name = 'PUNCHClient'):
    #     cls.__name__ = name

    # name_it(self)

    # @classmethod
    # def pre_search_hook(cls, *args, **kwargs):
    #     d = cls._get_match_dict(*args, **kwargs)
    #     # waverange = a.Wavelength(34*u.GHz, 17*u.GHz)
    #     # req_wave = d.get('Wavelength', waverange)
    #     # wmin = req_wave.min.to(u.GHz, equivalencies=u.spectral())
    #     # wmax = req_wave.max.to(u.GHz, equivalencies=u.spectral())
    #     # req_wave = a.Wavelength(wmin, wmax)
    #     # d['Wavelength'] = []
    #     # if 17*u.GHz in req_wave:
    #     #     d['Wavelength'].append('tca')
    #     # if 34*u.GHz in req_wave:
    #     #     d['Wavelength'].append('tcz')
    #     return cls.baseurl, cls.pattern, d


    @classmethod
    def register_values(cls):
        from sunpy.net import attrs as a

        adict = {attrs.Instrument: [('PUNCH', 'Polarimeter to UNify the Corona and Heliosphere')],
                attrs.Source: [('SwRI', 'Southwest Research Institute')],
                attrs.Provider: [('SwRI', 'Southwest Research Institute')],
                attrs.Level: [('3', 'level3')],
                PUNCHType: [('MPM', 'Magnetic Pseudo-Moments')],
                PUNCHLevel: [('level3', "Regular Level 3"),('level3_trefoil','The Trefoil Image')]
                }
        return adict
    
    @classmethod
    def _can_handle_query(cls, *query):
        required = {attrs.Instrument, attrs.Time}

        optional = {attrs.Level, attrs.Source, attrs.Provider, PUNCHType, PUNCHLevel}

        all_attrs = {type(x) for x in query}
        # [print(x, "\n\n") for x in all_attrs]
        return required.issubset(all_attrs) and all_attrs.issubset(required.union(optional))

    @classmethod
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
        from sunpy.net.scraper import Scraper
        from sunpy.time import TimeRange
        from sunpy.net.dataretriever.client import QueryResponse

        print("_______________________________________________________________")
        verb=True
        baseurl, pattern, matchdict = self.pre_search_hook(*args, **kwargs)
        # if verb: print("Args: ", args)
        # if verb: print("Baseurl: ", baseurl)
        # for ii, key in enumerate(matchdict.keys()):
        #     print(key," - ", matchdict[key])
        #     print(args[ii]," - ", args[ii].value)
        # if verb: print("Matchdict: ", matchdict)
        # if verb: print("Matchdict: ", matchdict["Instrument"][0])
        # if verb: print("SearchArg: ", args[0].value)
        tr = TimeRange(matchdict['Start Time'], matchdict['End Time'])

        filesmeta_list = []
        for key in matchdict:
            try:
                # if len(matchdict[key])>1:
                for value in matchdict[key]:
                    kwargs_dict = {key: value}
                    # print(kwargs_dict)
                    # print(key, value)
                    # print(baseurl)
                    this_baseurl = baseurl.replace("{{{}}}".format(key), value)
                    # print(this_baseurl)
                    # print("\n")
                    # import pdb; pdb.set_trace()
                    scraper = Scraper(this_baseurl, regex=True, **kwargs_dict)

                    filesmeta_list.extend(scraper._extract_files_meta(tr, extractor=pattern,
                                                matcher=matchdict))
                    # print(filesmeta_list)
                

            except TypeError:
                pass

        filesmeta = sorted(filesmeta_list, key=lambda k: k['url'])
        # print(filesmeta)
        # metalist = []
        # import pdb; pdb.set_trace()
        # for i in filesmeta:
        #     rowdict = self.post_search_hook(i, matchdict)
        #     metalist.append(rowdict)
        # print(metalist)

        return QueryResponse(filesmeta, client=self)


#run code if this is the main file
if __name__ == '__main__':
    from sunpy.net import Fido, attrs as a
    import os
    save_dir = os.getcwd()
    # print(save_dir)
    # r"%Y-%m-%dT%H:%M:%S"
    date_start = '2023-01-01T00:00:00'
    date_end = '2023-12-31T23:59:59'

    # print(a.Time(date_start, date_end))

    res = Fido.search(a.Instrument('punch'), a.Time(date_start, date_end)) #, PUNCHLevel('level3_trefoil'))
    print(res)
    # out = Fido.fetch(res, path=save_dir)
