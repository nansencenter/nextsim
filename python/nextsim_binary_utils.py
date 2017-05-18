import os,sys

def read_file_info(file_info):
      f  = open(file_info,'r')
      lines = f.readlines()
      f.close()

      variables       = []
      variable_types  = {} # list of variable types: 'f' or 'd' (single or double) or 'i' (integer)
      record_numbers  = {}

      for recno,lin in enumerate(lines):
         ss = lin.split()
         variables.append(ss[0])
         record_numbers.update({ss[0]:recno})
         if ss[1]=='double' or ss[1]=='8':
            variable_types.update({ss[0]:'d'})
         elif ss[1]=='single' or ss[1]=='float' or ss[1]=='4':
            variable_types.update({ss[0]:'f'})
         elif ss[1]=='int':
            variable_types.update({ss[0]:'i'}) # signed integer
         else:
            raise ValueError('unknown data type in '+file_info+': '+ss[1])

      return variables,record_numbers,variable_types

class nextsim_mesh_info:
   def __init__(self,mesh_file):

      self.mesh_file       = mesh_file
      self.mesh_file       = os.path.abspath(mesh_file)
      self.mesh_file_info  = self.mesh_file.replace('.bin','.dat')

      self.variables,self.record_numbers,self.variable_types  = read_file_info(self.mesh_file_info)

      return

class nextsim_binary_info:
   def __init__(self,data_file):

      # full or relative path to file
      self.data_file       = os.path.abspath(data_file)
      self.data_file_info  = self.data_file.replace('.bin','.dat')
      ss             = self.data_file.split('/')
      self.basename  = ss[-1].strip('.bin')
      self.dirname   = '/' #self.data_file.strip('/'+ss[-1])
      for sdir in ss[:-1]:
         self.dirname  += (sdir+'/')

      # field info
      self._read_field_info()
      self.variables_all   = self.variables

      # mesh info
      ss          = self.dirname+'/'+self.basename.replace('field','mesh')
      mesh_file   = ss +'.bin'
      self.mesh_info = nextsim_mesh_info(mesh_file)

      #TODO read nextsim.log (mppfile to initialise projection, original mesh_file)
      #TODO or pass them in as arguments?

      return

   def _read_field_info(self):

      import struct
      from datetime import datetime as DT,timedelta as TD

      # get the variables and their record no
      self.variables,self.record_numbers,self.variable_types  = read_file_info(self.data_file_info)

      # Time format
      fmt   = self.variable_types['Time']
      if fmt == 'f':
         size = 4
      elif fmt == 'd':
         size = 8

      f  = open(self.data_file,'rb')
      f.read(4)
      Time  = f.read(size)
      f.close()

      Time  = struct.unpack(fmt,Time)[0] #days
      
      refdate        = DT(1900,1,1)
      self.datetime  = refdate+TD(Time)
      self.datetimes = [self.datetime]

      return
