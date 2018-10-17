#!/usr/bin/python

import os
import time
import math
import ROOT as rt
from optparse import OptionParser # Command line parsing
import shutil
#from sets import Set
import random

def distance(a=(1,1),b=(2,2)):
   return math.sqrt((float(a[0])-float(b[0]))**2+(float(a[1])-float(b[1]))**2)

def empty_record():
   record = {}
   record["xSec"] = float(0.0)
   record["obsLimit"] = float(0.0)
   record["expLimit.m1sigma"] = float(0.0)
   record["expLimit.p1sigma"] = float(0.0)
   record["expLimit"] = float(0.0)
   return record

def main():
   usage = "usage: %prog options"
   version = "%prog."
   parser = OptionParser(usage=usage,version=version)
   parser.add_option("-i", "--inputdir", action="store", dest="inputdir", type="string", default=".", help="directory with mSugra files")
   parser.add_option("-n", "--noLog", action="store_true", dest="noLog", default=False)   
   parser.add_option("-g", "--grid", action="store", dest="grid", type="int", default = 25, help="interpolation grid")
   (options, args) = parser.parse_args()
   print "Input Directory: ", options.inputdir

   signal_key_list = []
   all_points = set([(0,0)])
   records = {}   

   print "Searching ", options.inputdir, " for mass points"

   for signal_name in os.listdir(options.inputdir):
      if not ("mSUGRA" in signal_name): continue
      if not (".dat" in signal_name): continue
      if signal_name.startswith('.'): continue
      if "tmp" in signal_name: continue
      record = {}
      record = empty_record()
      signal_file = open(options.inputdir+"/"+signal_name)

      mMother = -1
      mDaughter = -1

      for signal_line in signal_file:

         signal_line_split = signal_line.split()
         if "Mzero" in signal_line: 
            mMother = int(signal_line_split[-1])
            continue
         if "Mhalf" in signal_line: 
            mDaughter = int(signal_line_split[-1])
            continue
         if "Xsection" in signal_line:
            if options.noLog: 
               record["xSec"] = float(signal_line_split[-1])
            else:
               record["xSec"] = math.log(float(signal_line_split[-1]),10)
            continue

         if "error" in signal_line: continue
         if "1sigmaCoverage" in signal_line: continue
         if "2sigmaCoverage" in signal_line: continue
         if "m2sigma" in signal_line: continue
         if "p2sigma" in signal_line: continue

         if "observed" in signal_line in signal_line:
            obsRaw = float(signal_line_split[-1])
            if obsRaw <= 0: print mMother, mDaughter
            if options.noLog:
               record["obsLimit"] = float(signal_line_split[-1])
            else:
               record["obsLimit"] = math.log(float(signal_line_split[-1]),10)
            continue
         if "expected.m1sigma" in signal_line:
            if options.noLog: 
               record["expLimit.m1sigma"] = float(signal_line_split[-1])
            else:
               record["expLimit.m1sigma"] = math.log(float(signal_line_split[-1]),10)
            continue
         if "expected.p1sigma" in signal_line:
            if options.noLog:
               record["expLimit.p1sigma"] = float(signal_line_split[-1])
            else:
               record["expLimit.p1sigma"] = math.log(float(signal_line_split[-1]),10)
            continue
         if "expected" in signal_line:
            if options.noLog:
               record["expLimit"] = float(signal_line_split[-1])
            else:
               record["expLimit"] = math.log(float(signal_line_split[-1]),10)
            continue

#      print mMother+"_"+mDaughter
      if(mMother == -1 or mDaughter == -1):
         print signal_name, " has malformed data... skipping."
         continue

      if mDaughter == 1: mDaughter = 0
      signal_key_list.append((mMother,mDaughter))
      all_points.add((mMother,mDaughter))
      records[(mMother,mDaughter)] = record

      if mMother % 25 != 0: continue
      
      record = {}
      record = empty_record()

      grid = options.grid
#      grid = 25
#      grid = 50

      if mDaughter % grid == 0:
         dauStart = mDaughter
      else:
         dauStart = grid*((mDaughter//grid)+1)

      for i in range(0,dauStart,grid):
         if (mMother,i) not in all_points:
            all_points.add((mMother,i))
            records[(mMother,i)] = record
         if (mMother+grid,i) not in all_points:
            all_points.add((mMother+grid,i))
            records[(mMother+grid,i)] = record

   
   all_points.remove((0,0))

   given_points = set(signal_key_list)
 
   missing_points = all_points.difference(given_points)

   print "All points",len(all_points)
   print "Provided points",len(given_points)
   print "Missing points",len(missing_points)

   radius = 90

   neighbors = {}

   print "Determining nearest neighbors for missing points"

   radii = {}

   for miss in missing_points:
       close_neighbors = {}
       nearest_given = 60
       tmp_radius = radius
       given_count = 0
       for tmp in all_points:
           dist = distance(miss,tmp)      
           if dist < radius and dist > 1:
                close_neighbors[tmp] = 0
                if (tmp in given_points):
                   if (dist < nearest_given):
                      nearest_given = dist
                   given_count = given_count + 1

       while given_count > 4:
          given_count = 0
          tmp_radius = tmp_radius - (radius/20)
          remove_neighbors = []      

          for tmp in close_neighbors:
             if distance(miss,tmp) > tmp_radius: remove_neighbors.append(tmp)
             if tmp in given_points: given_count = given_count +1

          for tmp in remove_neighbors: del close_neighbors[tmp]

       radii[miss] = tmp_radius      

#       while given_count > 2:
       while given_count > 3:
          given_count = 0
          
#          farthest_distance = 0
          farthest_distance = 100000
          farthest_neighbor = miss

          for tmp in close_neighbors:
             if not tmp in given_points:
                continue

             given_count = given_count + 1             

#             if distance(miss,tmp)*random.uniform(.95,1.05) > farthest_distance:
#             if distance((0,0),tmp)*random.uniform(.99,1.01) > farthest_distance:
#Let's try being agressive in the opposite direction
             if distance((0,0),tmp)*random.uniform(.99,1.01) < farthest_distance:
                farthest_distance = distance((0,0),tmp)*random.uniform(.99,1.01)
#                farthest_distance = distance(miss,tmp)*random.uniform(.95,1.05)
                farthest_neighbor = tmp

          del close_neighbors[farthest_neighbor]
          given_count = given_count - 1

       neighbors[miss] = close_neighbors

   #print radii
#Let's start with some reasonably good values
   print "Seeding missing points with reasonable values"
   for miss in missing_points:
      best_distance = 10000
      tmp_point = (-1,-1)
      for tmp in given_points:
         if distance(miss,tmp)*random.uniform(.95,1.05) < best_distance:
            best_distance = distance(miss,tmp)*random.uniform(.95,1.05)
            tmp_point = tmp
 
      records[miss]["xSec"] = records[tmp_point]["xSec"]*random.uniform(.95,1.05)
      records[miss]["obsLimit"] = records[tmp_point]["obsLimit"]*random.uniform(.95,1.05)
      records[miss]["expLimit.m1sigma"] = records[tmp_point]["expLimit.m1sigma"]*random.uniform(.95,1.05)
      records[miss]["expLimit.p1sigma"] = records[tmp_point]["expLimit.p1sigma"]*random.uniform(.95,1.05)
      records[miss]["expLimit"] = records[tmp_point]["expLimit"]*random.uniform(.95,1.05) 

   print "Determining neighbor weights for missing points"

   count = 0
   for miss in missing_points:
      count = count + 1
      if count % 25 == 0: print count,"points processed of",len(missing_points)
      sum = 0
#      step = float(float(radius) * 1.05 / 15)
#      for i in range(-15,16):
#         for j in range(-15,16):
#            x = float(miss[0])+float(i)*step
#            y = float(miss[1])+float(j)*step
#            print miss, (x,y)
#            if x==miss[0] and y==miss[1]: continue
#            if distance((x,y),miss) > radius: continue
      for i in range(0,1000):
         r = random.uniform(0,radii[miss])
         theta = random.uniform(-math.pi,math.pi)
         x = miss[0]+r*math.cos(theta)
         y = miss[1]+r*math.sin(theta)
         if x != miss[0] and y != miss[1]:
            best_distance = 1000
            best_neighbor = (-1,-1)
            for neighbor in neighbors[miss]:
               if (distance((x,y),neighbor)*random.uniform(.95,1.05)) < best_distance:
                  best_distance = (distance((x,y),neighbor)*random.uniform(.95,1.05))
                  best_neighbor = neighbor
            if best_neighbor in given_points:
               contrib = 10.0
            else:
               contrib = 1.0
            neighbors[miss][best_neighbor] = float(neighbors[miss][best_neighbor]) + contrib
            sum = sum + contrib

      for neighbor in neighbors[miss]: neighbors[miss][neighbor] = float(neighbors[miss][neighbor]) / sum
   print "Neighbor weights have been calculated"

   print "Diffusing..."
   error = 1
   count = 0
   while error > .0000001 or count < 250:
      count = count + 1
      if count % 25 == 0 : print "   diffuse step",count
      new_records = {}
      error = 0
      for miss in missing_points:
         record = empty_record()
         for neighbor in neighbors[miss]:
            record["xSec"] = record["xSec"] + (records[neighbor]["xSec"] * neighbors[miss][neighbor])
            record["obsLimit"] = record["obsLimit"] + (records[neighbor]["obsLimit"] * neighbors[miss][neighbor])
            record["expLimit.m1sigma"] = record["expLimit.m1sigma"] + (records[neighbor]["expLimit.m1sigma"] * neighbors[miss][neighbor])
            record["expLimit.p1sigma"] = record["expLimit.p1sigma"] + (records[neighbor]["expLimit.p1sigma"] * neighbors[miss][neighbor])
            record["expLimit"] = record["expLimit"] + (records[neighbor]["expLimit"] * neighbors[miss][neighbor])
         new_records[miss] = record
#         if (records[miss]["expLimit"] + new_records[miss]["expLimit"]) <= 0: continue
         tmp = 2 * abs(records[miss]["expLimit"] - new_records[miss]["expLimit"]) / (records[miss]["expLimit"] + new_records[miss]["expLimit"])
         if tmp > error: error = tmp 
      
      for miss in missing_points:
         records[miss] = new_records[miss]
      
      #print "Max error:",error

   print "Diffusing complete"

   print "Saving Interpolated Files"
   for miss in missing_points:
      momMass = str(miss[0])
      if miss[1] == 0:
         dauMass = "1"
      else:
         dauMass = str(miss[1])

      format_file_name = "mSUGRA_" + momMass + "_" + dauMass + "_10_0_1.dat"

      format_file = open(options.inputdir+"/"+format_file_name, "w")

      if options.noLog:
         strXsec     = str(records[miss]["xSec"])
         strLimObs   = str(records[miss]["obsLimit"])
         strLimExpM  = str(records[miss]["expLimit.m1sigma"])
         strLimExp   = str(records[miss]["expLimit"])
         strLimExpP  = str(records[miss]["expLimit.p1sigma"])

      else:
         strXsec     = str(math.pow(10,records[miss]["xSec"]))
         strLimObs   = str(math.pow(10,records[miss]["obsLimit"]))
         strLimExpM  = str(math.pow(10,records[miss]["expLimit.m1sigma"]))
         strLimExp   = str(math.pow(10,records[miss]["expLimit"]))
         strLimExpP  = str(math.pow(10,records[miss]["expLimit.p1sigma"]))


      formatLine = """Azero = 0
Mu = 1
Mzero = """ + momMass + """
Mhalf = """ + dauMass + """
Xsection = """ + strXsec + """
limit.cls.observed = """ + strLimObs + """
limit.cls.observed.error = 0
limit.cls.expected.m1sigma = """ + strLimExpM + """
limit.cls.expected = """ + strLimExp + """
limit.cls.expected.p1sigma = """ + strLimExpP + """
limit.cls.expected.1sigmaCoverage = -1
limit.cls.expected.2sigmaCoverage = -1
limit.cls.expected.m2sigma = -1
limit.cls.expected.p2sigma = -1

"""

      format_file.write(formatLine)
      format_file.close()

   exit()

if __name__ == "__main__":
   main()

