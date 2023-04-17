# -*- coding: utf-8 -*-
'''
PyCORN - script to extract data from .res (results) files generated
by UNICORN Chromatography software supplied with ÄKTA Systems
(c)2014-2016 - Yasar L. Ahmed
v0.18b
'''

from __future__ import print_function
from collections import OrderedDict
from zipfile import ZipFile
from zipfile import is_zipfile
import xml.etree.ElementTree as ET
import struct
import codecs
import os
import io
import pandas as pd
from copy import deepcopy

class pc_res3(OrderedDict):
    """A class for holding the PyCORN/RESv3 data.
    A subclass of `dict`, with the form `data_name`: `data`.
    """

    # first, some magic numbers
    RES_magic_id = b'\x11\x47\x11\x47\x18\x00\x00\x00\xB0\x02\x00\x00\x20\x6C\x03\x00'
    CNotes_id = b'\x00\x00\x01\x00\x02\x00\x03\x22'
    Methods_id = b'\x00\x00\x01\x00\x02\x00\x01\x02'
    Logbook_id = b'\x00\x00\x01\x00\x04\x00\x48\x04'
    Logbook_id2 = b'\x00\x00\x01\x00\x04\x00\x49\x04'
    SensData_id = b'\x00\x00\x01\x00\x04\x00\x01\x14'
    SensData_id2 = b'\x00\x00\x01\x00\x04\x00\x02\x14'
    Fractions_id = b'\x00\x00\x01\x00\x04\x00\x44\x04'
    Fractions_id2 = b'\x00\x00\x01\x00\x04\x00\x45\x04'
    Inject_id = b'\x00\x00\x01\x00\x04\x00\x46\x04'
    Inject_id2 = b'\x00\x00\x01\x00\x04\x00\x47\x04'
    LogBook_id = b'\x00\x00\x01\x00\x02\x00\x01\x13'  # capital B!

    def __init__(self, file_name, reduce=1, inj_sel=-1):
        OrderedDict.__init__(self)
        self.file_name = file_name
        self.reduce = reduce
        self.injection_points = None
        self.inj_sel = inj_sel
        self.inject_vol = None
        self.header_read = False
        self.run_name = ''

        with open(self.file_name, 'rb') as f:
            self.raw_data = f.read()

    def input_check(self, show=False):
        '''
        Checks if input file is a supported res file
        x = magic number, y = version string, z = file size/EOF

        Returns True or False
        '''
        if show: print((" ---- \n Input file: {0}").format(self.file_name))

        x = self.raw_data.find(self.RES_magic_id, 0, 16)
        y = self.raw_data.find(b'UNICORN 3.10', 16, 36)
        z = struct.unpack("i", self.raw_data[16:20])

        if (x, y) == (0, 24):
            if show: print(" Input is a UNICORN 3.10 file!")
            x, y = (0, 0)
        else:
            if show: print(" Input is not a UNICORN 3.10 file!")
            x, y = (1, 1)

        if z[0] == os.path.getsize(self.file_name):
            if show: print(" File size check - OK")
            z = 0
        else:
            if show: print(" File size mismatch - file corrupted?")
            z = 1
        if (x, y, z) != (0, 0, 0):
            if show: print("\n File not supported - stop!")
            return False
        else:
            if show: print("\n Alles safe! - Go go go!")
            return True

    def readheader(self):
        '''
        Extracts all the entries/declarations in the header (starts at position 686)
        '''

        # we only need to do this once
        if self.header_read: return
        self.header_read = True

        fread = self.raw_data
        header_end = fread.find(self.LogBook_id) + 342
        for i in range(686, header_end, 344):
            decl = struct.unpack("8s296s4i", fread[i:i + 320])
            full_label = codecs.decode(decl[1], 'iso8859-1').rstrip("\x00")
            if full_label.find(':') == -1:
                r_name = ''
                d_name = full_label
            else:
                r_name = full_label[:full_label.find(':')]
                d_name = full_label[full_label.find('_') + 1:]
            x = dict(magic_id=decl[0],
                     run_name=r_name,
                     data_name=d_name,
                     d_size=decl[2],
                     off_next=decl[3],
                     adresse=decl[4],
                     off_data=decl[5],
                     d_start=decl[4] + decl[5],
                     d_end=decl[4] + decl[2])
            name = x['data_name']
            dat = self.get(name, dict())
            dat.update(x)
            self[name] = dat

    def showheader(self, full=True):
        '''
        Prints content of header
        '''
        print((" ---- \n Header of {0}: \n").format(self.file_name))
        if full:
            print("  MAGIC_ID, ENTRY_NAME, BLOCK_SIZE, OFFSET_TO_NEXT, ADRESSE, OFFSET_TO_DATA")
        else:
            print(" ENTRY_NAME, BLOCK_SIZE, OFFSET_TO_NEXT, ADRESSE, OFFSET_TO_DATA")
        num_blocks = len(self.items())
        for i in range(num_blocks):
            dtp = (list(self.items()))[i][1]
            if full:
                print(" ", dtp['magic_id'], dtp['data_name'], dtp['d_size'], dtp['off_next'], dtp['adresse'],
                      dtp['off_data'])
            else:
                print(" ", dtp['data_name'], dtp['d_size'], dtp['off_next'], dtp['adresse'], dtp['off_data'])

    def get_user(self):
        '''
        Show stored user name
        '''
        fread = self.raw_data[:512]
        u = struct.unpack("40s", fread[118:158])
        dec_u = codecs.decode(u[0], 'iso8859-1').rstrip("\x00")
        return dec_u

    def dataextractor(self, dat, show=False):
        '''
        Identify data type by comparing magic id, then run appropriate
        function to extract data, update orig. dict to include new data
        '''
        meta1 = [
            self.Logbook_id, self.Logbook_id2,
            self.Inject_id, self.Inject_id2,
            self.Fractions_id, self.Fractions_id2]
        meta2 = [self.CNotes_id, self.Methods_id]
        sensor = [self.SensData_id, self.SensData_id2]
        if dat['d_size'] == 0:
            pass
        elif dat['magic_id'] in meta1:
            dat.update(data=self.meta1_read(dat, show=show), data_type= 'annotation')
            return dat
        elif dat['magic_id'] in meta2:
            dat.update(data=self.meta2_read(dat, show=show), data_type= 'meta')
            return dat
        elif dat['magic_id'] in sensor:
            values, unit = self.sensor_read(dat, show=show)
            dat.update(data=values, unit=unit, data_type= 'curve')
            return dat

    def meta1_read(self, dat, show=False, do_it_for_inj_det=False):
        '''
        Extracts meta-data/type1, Logbook, fractions and Inject marks
        for a specific datum
        '''
        if show:
            print((" Reading: {0}").format(dat['data_name']))
        final_data = []
        inj_vol_to_subtract = self.inject_vol
        if do_it_for_inj_det:
            inj_vol_to_subtract = 0.0     
        for i in range(dat['d_start'], dat['d_end'], 180):
            dp = struct.unpack("dd158s", self.raw_data[i:i + 174])
            # acc_time = dp[0] # not used atm
            acc_volume = round(dp[1] - inj_vol_to_subtract, 4)
            label = (codecs.decode(dp[2], 'iso8859-1')).rstrip('\x00')
            merged_data = acc_volume, label
            final_data.append(merged_data)
        return (final_data)

    def meta2_read(self, dat, show=False):
        '''
        Extracts meta-data/type2, Method/Program used in the run
        '''
        if show: print((" Reading: {0}").format(dat['data_name']))
        start, size = dat['d_start'], dat['d_size']
        tmp_data = self.raw_data[start:start + size]
        size = tmp_data.rfind(b'\n')  # declared block-size in header is always off
        # by a few bytes, hence it is redetermined here
        if show and size != len(tmp_data):
            print('meta2: reevaluated size {} -> {}'.format(size, len(tmp_data)))

        raw_data = codecs.decode(self.raw_data[start:start + size], 'iso8859-1')
        if '\r' in raw_data:
            data = raw_data
        else:
            data = raw_data.replace('\n', '\r\n')
        return data

    def sensor_read(self, dat, show=False):
        '''
        extracts sensor/run-data and applies correct division
        '''
        final_data = []
        if "UV" in dat['data_name'] or "Cond" == dat['data_name'] or "Flow" == dat['data_name']:
            sensor_div = 1000.0
        elif "Pressure" in dat['data_name']:
            sensor_div = 100.0
        else:
            sensor_div = 10.0
        if show: print((" Reading: {0}").format(dat['data_name']))

        fread = self.raw_data
        for i in range(dat['adresse'] + 207, dat['adresse'] + 222, 15):
            s_unit = struct.unpack("15s", fread[i:i + 15])
            s_unit_dec = (codecs.decode(s_unit[0], 'iso8859-1')).rstrip('\x00')
            # FIX: in some files the unit for temperature reads 'C' instead of '°C' 
            if s_unit_dec == 'C':
                s_unit_dec = u'°C'
        for i in range(dat['d_start'], dat['d_end'], 8):
            sread = struct.unpack("ii", fread[i:i + 8])
            data = round((sread[0] / 100.0) - self.inject_vol, 4), sread[1] / sensor_div
            final_data.append(data)
        return (final_data[0::self.reduce], s_unit_dec)

    def inject_det(self, show=False):
        '''
        Finds injection points - required for adjusting retention volume
        '''
        inject_ids = [self.Inject_id, self.Inject_id2]
        injections = []
        if self.injection_points == None:
            self.injection_points = [0.0]
            for i in self.values():
                if i['magic_id'] in inject_ids:
                    injections = self.meta1_read(i, show=show, do_it_for_inj_det=True)
        for i in injections:
            if i[0] != 0.0:
                self.injection_points.append(i[0])
        if show:
            print(" ---- \n Injection points: \n # \t ml")
            for x, y in enumerate(self.injection_points):
                print((" {0} \t {1}").format(x, y))


    def load(self, show=False):
        '''
        extract all data and store in list
        '''
        self.readheader()
        self.run_name = self['Logbook']['run_name']
        self.inject_det()
        try:
            self.inject_vol = self.injection_points[self.inj_sel]
        except IndexError:
            print("\n WARNING - Injection point does not exist! Selected default.\n")
            self.inject_vol = self.injection_points[-1]
        for name, dat in list(self.items()):
            dat = self.dataextractor(dat, show=show)
            if dat is not None:
                self[name] = dat
            else:
                # TODO: Maybe we should keep this around?
                del self[name]
                
class pc_uni6(OrderedDict):
    '''
    A class for holding the pycorn/RESv6 data
    A subclass of `dict`, with the form `data_name`: `data`.
    '''
    # for manual zip-detection
    zip_magic_start = b'\x50\x4B\x03\x04\x2D\x00\x00\x00\x08'
    zip_magic_end = b'\x50\x4B\x05\x06\x00\x00\x00\x00'
    
    # hack to get pycorn-bin to move on
    SensData_id = 0
    SensData_id2 = 0
    Fractions_id = 0
    Fractions_id2 = 0
    
    def __init__(self, inp_file, r_name='0p1l2a3c4e5h6o7l8d9e0r', inj_vol=0.0):
        OrderedDict.__init__(self)
        if r_name == '0p1l2a3c4e5h6o7l8d9e0r':
            self.run_name = inp_file
        else:
            self.run_name = r_name
        self.file_name = inp_file
        self.load()
        self.xml_parse()
        self.clean_up()
        if inj_vol == 0.0 and 'Injection' in self and 'data' in self['Injection']:
            self.inject_vol = self['Injection']['data'][0][0]
        else:
            self.inject_vol = inj_vol
        #self.substract_inj_vol()

    
    def load(self, show=False):
        '''
        zip-files inside the zip-bundle are replaced by dicts, again with dicts with filename:content
        Chrom.#_#_True (=zip-files) files are unpacked from binary to floats by unpacker()
        To access x/y-value of Chrom.1_2:
        udata = pc_uni6("mybundle.zip")
        udata.load()
        x = udata['Chrom.1_2_True']['CoordinateData.Volumes']
        y = udata['Chrom.1_2_True']['CoordinateData.Amplitudes']
        '''
        with open(self.file_name, 'rb') as f:
            input_zip = ZipFile(f)
            zip_data = self.zip2dict(input_zip)
            self.update(zip_data)
            proc_yes = []
            proc_no = []
            for i in self.keys():
                tmp_raw = io.BytesIO(input_zip.read(i))
                f_header = tmp_raw.read(9)
                # tmp_raw.seek(0)
                # the following if block is to fix the non-standard zip files
                # by stripping out all the null-bytes at the end
                # see https://bugs.python.org/issue24621
                if f_header == self.zip_magic_start:
                    proper_zip = tmp_raw.getvalue()
                    f_end = proper_zip.rindex(self.zip_magic_end) + 22
                    tmp_raw = io.BytesIO(proper_zip[0:f_end])
                if is_zipfile(tmp_raw):
                    tmp_zip = ZipFile(tmp_raw)
                    x = {i:self.zip2dict(tmp_zip)}
                    self.update(x)
                    proc_yes.append(i)
                else:
                    pass
                    proc_no.append(i)
            if show:
                print("Loaded " + self.file_name + " into memory")
                print("\n-Supported-")
                for i in proc_yes:
                    print(" " + i)
                print("\n-Not supported-")
                for i in proc_no:
                    print(" " + i)
        # filter out data we dont deal with atm
        to_process = []
        for i in self.keys():
            if "Chrom" in i and not "Xml" in i:
                to_process.append(i)
        if show:
            print("\nFiles to process:")
            for i in to_process:
                print(" " + i)
        for i in to_process:
            for n in self[i].keys():
                if "DataType" in n:
                    a = self[i][n]
                    b = a.decode('utf-8')
                    x = b.strip("\r\n")
                else:
                    x = self.unpacker(self[i][n])
                tmp_dict = {n:x}
                self[i].update(tmp_dict)
        if show:
            print("Finished decoding x/y-data!")

    @staticmethod
    def zip2dict(inp):
        '''
        input = zip object
        output = dict with filename:file-object pairs
        '''
        mydict = {}
        for i in inp.NameToInfo:
            tmp_dict = {i:inp.read(i)}
            mydict.update(tmp_dict)
        return(mydict)
    
    @staticmethod
    def unpacker(inp):
        '''
        input = data block
        output = list of values
        '''
        read_size = len(inp) - 48
        values = []
        for i in range(47, read_size, 4):
            x = struct.unpack("<f", inp[i:i+4])
            x = x[0]
            values.append(x)
        return(values)
   
    def xml_parse(self,show=False):
        '''
        parses parts of the Chrom.1.Xml and creates a res3-like dict
        '''
        if 'Chrom.1.Xml' in self.keys():
            tree = ET.fromstring(self['Chrom.1.Xml'])
        elif 'Compare.Xml' in self.keys():
            tree = ET.fromstring(self['Compare.Xml'])
        mc = tree.find('Curves')
        me = tree.find('EventCurves')
        print(tree.tag)
        print(tree.attrib)
        event_dict = {}
        for i in range(len(me)):
            magic_id = self.SensData_id
            e_type = me[i].attrib['EventCurveType']
            e_name = me[i].find('Name').text
            if e_name == 'Fraction':
                e_name = 'Fractions' # another hack for pycorn-bin
            e_orig = me[i].find('IsOriginalData').text
            e_list = me[i].find('Events')
            e_data = []
            for e in range(len(e_list)):
                e_vol = float(e_list[e].find('EventVolume').text)
                e_txt = e_list[e].find('EventText').text
                e_data.append((e_vol,e_txt))
            if e_orig == "true":
                # substract injection volume from x data (volume data)
                injcurvenumber = mc[i].find('./InjectionEventCurve/CurveNumber')
                if injcurvenumber is not None:
                    injcurvenumber = injcurvenumber.text
                    injvolume = me.find('./EventCurve/CurveNumber/[.="' + injcurvenumber +  '"]..Events/*[@EventType="Injection"]/EventVolume').text
                    injvolume = float(injvolume)
                    e_data = [tuple([e[0] - injvolume, ] + list(e[1:])) for e in e_data]
                x = {'run_name':self.run_name, 'data': e_data, 'data_name':e_name, 'magic_id':magic_id}
                event_dict.update({e_name:x})
        self.update(event_dict)
        chrom_dict = {}
        for i in range(len(mc)):
            d_type = mc[i].attrib['CurveDataType']
            d_name = mc[i].find('Name').text
            d_fname = mc[i].find('CurvePoints')[0][1].text # equals mc[i].find('CurvePoints/CurvePoint/BinaryCurvePointsFileName').text
            d_unit = mc[i].find('AmplitudeUnit').text
            magic_id = self.SensData_id
            if 'CoordinateData.Amplitudes' in self[d_fname].keys():
                y_dat = self[d_fname]['CoordinateData.Amplitudes']
                x_dat = None

                if 'CoordinateData.Volumes' in self[d_fname].keys():
                    x_dat = self[d_fname]['CoordinateData.Volumes']
                elif mc[i].find('DerivedCurveOrigin') is not None:
                    # curves that are derived from other curves dont have x data
                    # primcurvenumber = mc[i].find('DerivedCurveOrigin/FirstCurve/CurveNumber').text
                    # primcurvename = mc.find('./Curve/CurveNumber/[.="' + primcurvenumber +  '"]..CurvePoints/CurvePoint/BinaryCurvePointsFileName').text
                    # primcurvename = mc.find('./Curve/Name/[.="System flow"]..CurvePoints/CurvePoint/BinaryCurvePointsFileName').text
                    startpoint = mc[i].find('DistanceToStartPoint').text
                    startpoint = float(startpoint)
                    pointdistance = mc[i].find('DistanceBetweenPoints').text
                    pointdistance = float(pointdistance)
                    x_dat = []
                    for j in range(len(y_dat)):
                        x_dat.append(startpoint + j * pointdistance)
                else:
                    pass
                    # don't deal with data that does not make sense atm
                    # orig2.zip contains UV-blocks that are (edited) copies of
                    # original UV-trace but they dont have the volume data

                if x_dat is not None:
                    # substract injection volume from x_dat (volume data)
                    injcurvenumber = mc[i].find('./InjectionEventCurve/CurveNumber')
                    if injcurvenumber is not None:
                        injcurvenumber = injcurvenumber.text
                        injvolume = me.find('./EventCurve/CurveNumber/[.="' + injcurvenumber +  '"]..Events/*[@EventType="Injection"]/EventVolume').text
                        injvolume = float(injvolume)
                        for x in range(len(x_dat)):
                            x_dat[x] = x_dat[x] - injvolume
                    zdata = list(zip(x_dat,y_dat))
                    # hack to prevent pycorn-bin from picking this up
                    if d_name == "UV cell path length":
                        d_name = "xUV cell path length"
                    # data type for flow rate
                    if d_unit == "ml/min":
                        d_type = "Flow rate"
                    x = {'run_name':self.run_name, 'data': zdata, 'unit': d_unit, 'data_name':d_name, 'data_type':d_type, 'magic_id':magic_id}
                    chrom_dict.update({d_name:x})
            if show:
                print("---")
                print(d_type)
                print(d_name)
                print(d_fname)
                print(d_unit)
        self.update(chrom_dict)
    
    def clean_up(self):
        '''
        deletes everything and just keeps relevant run-date
        resulting dict is more like res3
        '''
        manifest = ET.fromstring(self['Manifest.xml'])
        for i in range(len(manifest)):
            file_name = manifest[i][0].text
            self.pop(file_name)
        self.pop('Manifest.xml')

    def xy_data(self, key):
        '''
        Takes a data block and returns two lists with x- and y-data
        '''
        if 'data' in self[key]:
            inp = self[key]['data']
            x_data = [x[0] for x in inp]
            y_data = [x[1] for x in inp]
            return x_data, y_data

    def write_csv(self, fname):
        '''
        writes sensor/run-data to csv-files
        '''
        for i in self.keys():
            print("Writing: " + self[i]['data_name'])
            outfile_base = fname[:-4] + "_" + self.run_name + "_" + self[i]['data_name']
            type = self[i]['data_type']
            if type == 'meta':
                data = self[i]['data']
                data_to_write = data.encode('utf-8')
                ext = '.txt'
                sep = '\t'
                with open(outfile_base + ext, 'wb') as fout:
                    fout.write(data_to_write)
            else:
                x_dat,y_dat = self.xy_data(i)
                ext = '.csv'
                sep = ','
                with open(outfile_base + ext, 'wb') as fout:
                    for x,y in zip(x_dat,y_dat):
                        dp = str(x) + sep + str(y) + str('\r\n')
                        data_to_write = dp.encode('utf-8')
                        fout.write(data_to_write)

    def write_xls(self, fname):
        '''
        Input = pycorn object
        output = xlsx file
        '''
        import xlsxwriter
        xls_filename = fname[:-4] + "_" + self.run_name + ".xlsx"
        workbook = xlsxwriter.Workbook(xls_filename)
        worksheet = workbook.add_worksheet()
        writable_blocks = [self.Fractions_id, self.Fractions_id2, self.SensData_id, self.SensData_id2]
        d_list = []
        for i in self.keys():
            if self[i]['magic_id'] in writable_blocks:
                d_list.append(i)
        for i in d_list:
            dat = self[i]['data']
            try:
                unit = self[i]['unit']
            except:
                KeyError
                unit = 'Fraction'
            header1 = (self[i]['data_name'], '')
            header2 = ('ml', unit)
            dat.insert(0, header1)
            dat.insert(1, header2)
            row = 0
            col = d_list.index(i) *2
            print("Writing: " + i)
            for x_val, y_val in (dat):
                worksheet.write(row, col, x_val)
                worksheet.write(row, col + 1, y_val)
                row += 1
        workbook.close()
        print("Data written to: " + xls_filename)

    def to_pandas(self):
        '''
        converts the dict to a pandas dataframe
        '''
        df = pd.DataFrame.from_dict(self, orient='index')
        df = df.explode('data')
        df[['volume','value']] = pd.DataFrame(df['data'].tolist(), index=df.index)
        df = df.drop(columns=['data'])
        df['volume'] = pd.to_numeric(df['volume'])

        return(df)
    
    def volume_to_time(self, key):
        '''
        Returns a copy of self[key] with volume data converted to time data based on system flow
        '''
        if 'data' in self[key]:
            v_x_dat, y_dat = self.xy_data('System flow')
            t_x_dat = [0,]
            # calculate time points based on system flow
            for i in range(len(v_x_dat)-1):
                if y_dat[i] == 0:
                    dt_v = 0
                else:
                    dt_v = (v_x_dat[i+1] - v_x_dat[i])/y_dat[i]
                if i == 0:
                    t_v = dt_v
                else:
                    t_v = t_x_dat[i] + dt_v
                t_x_dat.append(t_v)

            # keep time points before injection negative
            inj_index = v_x_dat.index(min(v_x_dat, key=lambda x:abs(x)))
            t_x_dat = [x - t_x_dat[inj_index] for x in t_x_dat]

            # copy data block for conversion
            data = deepcopy(self[key])

            # convert volume data to time data of requested data block
            p_t_x_dat = []
            for p in data['data']:
                closest_1_key = 1
                closest_2_key = 0
                # find the two closest volume points
                for i in range(2, len(t_x_dat)):
                    if abs(p[0] - v_x_dat[i]) <= abs(p[0] - v_x_dat[closest_1_key]):
                        closest_2_key = closest_1_key
                        closest_1_key = i
                    elif abs(p[0] - v_x_dat[i+1]) <= abs(p[0] - v_x_dat[closest_2_key]):
                        closest_2_key = i+1
                    else:
                        break

                # calculate time point based on linear interpolation
                m = (t_x_dat[closest_2_key] - t_x_dat[closest_1_key])/(v_x_dat[closest_2_key] - v_x_dat[closest_1_key])
                b = t_x_dat[closest_1_key] - m*v_x_dat[closest_1_key]
                t = m*p[0] + b
                p_t_x_dat.append(t)
            
            # replace volume data with time data
            data['data'] = [tuple([p_t_x_dat[i], ] + list(data['data'][i][1:])) for i in range(len(data['data']))]
            return data