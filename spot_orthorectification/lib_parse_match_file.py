#!/usr/bin/env python
# __BEGIN_LICENSE__
#  Copyright (c) 2009-2013, United States Government as represented by the
#  Administrator of the National Aeronautics and Space Administration. All
#  rights reserved.
#
#  The NGT platform is licensed under the Apache License, Version 2.0 (the
#  "License"); you may not use this file except in compliance with the
#  License. You may obtain a copy of the License at
#  http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
# __END_LICENSE__

# Thanks to Amaury Dehecq for contributing this tool.
# Bugfix by PicoJr.
# Transformed into a library by Enrico Mattea (2023/01/25).

import os, struct
import numpy as np

def read_ip_record(mf):
    """
    Read one IP record from the binary match file.
    Information comtained are x, y, xi, yi, orientation, scale, interest, polarity, octave, scale_lvl, desc 
    Input: - mf, file handle to the in put binary file (in 'rb' mode)
    Output: - iprec, array containing the IP record
    """
    x, y = np.frombuffer(mf.read(8), dtype=np.float32)
    xi, yi = np.frombuffer(mf.read(8), dtype=np.int32)
    orientation, scale, interest = np.frombuffer(mf.read(12), dtype=np.float32)
    polarity, = np.frombuffer(mf.read(1), dtype=np.int8)  # or np.bool?
    octave, scale_lvl = np.frombuffer(mf.read(8), dtype=np.uint32)
    ndesc, = np.frombuffer(mf.read(8), dtype=np.uint64)
    desc = np.frombuffer(mf.read(int(ndesc * 4)), dtype=np.float32)
    iprec = [x, y, xi, yi, orientation, scale, interest, polarity, octave, scale_lvl, ndesc]
    iprec.extend(desc)
    return iprec


def read_match_file(match_file):
    """
    Read a full binary match file. First two 8-bits contain the number of IPs in each image. Then contains the record for each IP, image1 first, then image2.
    Input: 
    - match_file: str, path to the match file
    Outputs:
    - two arrays, containing the IP records for image1 and image2.
    """

    # Open binary file in read mode
    print("Reading: " + match_file)
    mf = open(match_file,'rb')

    # Read record length
    size1 = np.frombuffer(mf.read(8), dtype=np.uint64)[0]
    size2 = np.frombuffer(mf.read(8), dtype=np.uint64)[0]

    # Read record for each image
    im1_ip = [read_ip_record(mf) for i in range(size1)]
    im2_ip = [read_ip_record(mf) for i in range(size2)]

    # Close file
    mf.close()
    
    return im1_ip, im2_ip



def process(infile, outfile):

    # Read match file
    im1_ip, im2_ip = read_match_file(infile)


    # Extract only interest point positions, combine point-wise.
    l1=[sublist[:2] for sublist in im1_ip]
    l2=[sublist[:2] for sublist in im2_ip]

    ip_out = [l1[i]+l2[i] for i in range(len(l1))]

    # Save to text file
    print("Writing: " + outfile)

    try:
        os.makedirs(os.path.dirname(outfile))
    except OSError:
        pass
        
    with open(outfile, 'w') as out:

        # Write IPs
        for i in range(len(ip_out)):
            iprec = ip_out[i]
            iprec_str = [str(a) for a in iprec]
            out.write('\t'.join(iprec_str) + '\n')
