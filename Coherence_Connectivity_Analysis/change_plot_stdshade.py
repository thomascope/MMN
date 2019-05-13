#!/usr/bin/python
"""
Created on Fri Apr 26 16:20:59 2019
A script for swapping 'plot' commands to stdshade commands, assuming all 
averaging across the specified dimension (for MMN granger and icoh, the sixth, subject number)
@author: tc02
"""

import re
pattern = re.compile("plot\(foi,squeeze\(.*,.*\)")
oldfunc = 'plot';
newfunc = 'stdshade_TEC';
#pattern = "*plot*";

with open('scratch_compare.m') as file_for_swapping, open('scratch_compare_substituted.m','w') as newfile:
    for line in file_for_swapping:
        found_string = pattern.search(line);
        if found_string and 'random' not in line:
            #print "Found this line: ", line
            this_string = found_string.group()
            this_string = this_string.replace(oldfunc,newfunc)
            this_split_string = this_string.split('(',1)            
            this_split_string_after = this_split_string[1].split(',')
            try:            
                if len(this_split_string_after) == 12:   
                    if 'mean' in this_split_string_after[1]:
                        new_func = this_split_string[0] + '(' + this_split_string_after[1].replace('mean(','',1) + ',' + ','.join(this_split_string_after[2:8]) + ')\',0.2,' +  this_split_string_after[9] + ',' + this_split_string_after[0] + ',1,1)'
                    elif 'median' in this_split_string_after[1]:
                        new_func = this_split_string[0] + '(' + this_split_string_after[1].replace('median(','',1) + ',' + ','.join(this_split_string_after[2:8]) + ')\',0.2,' +  this_split_string_after[9] + ',' + this_split_string_after[0] + ',1,1)'
                elif len(this_split_string_after) == 10: 
                    if 'mean' in this_split_string_after[1]:
                        new_func = this_split_string[0] + '(' + this_split_string_after[1].replace('mean(','',1) + ',' + ','.join(this_split_string_after[2:6]) + ')\',0.2,' +  this_split_string_after[7] + ',' + this_split_string_after[0] + ',1,1)'
                    elif 'median' in this_split_string_after[1]:
                        new_func = this_split_string[0] + '(' + this_split_string_after[1].replace('median(','',1) + ',' + ','.join(this_split_string_after[2:6]) + ')\',0.2,' +  this_split_string_after[7] + ',' + this_split_string_after[0] + ',1,1)'
                line = line.replace(found_string.group(),new_func)
            except:
                #print line
                pass
        newfile.write(line)
            