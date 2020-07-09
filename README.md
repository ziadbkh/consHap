
# consHap 
**consHap (v0.1)**: Majority voting approach to construct consensus haplotype estimator from multiple phased haplotypes.

# Licence
Copyright **2020 Ziad Al Bkhetan**

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at

```
http://www.apache.org/licenses/LICENSE-2.0
```

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.

# Prerequisites

 1. Python 3.5.2 or later.


# Commands
## Help

    python3 consHap.py --help

## Consensus haplotype construction 
    python3 consHap.py -CONS -O path_output_file_prefix -F list_of_input_files 

You can use the option **-GZ** for gzipped input and output files.    
You can use the option **-RI** to impute missing SNPs as the reference file (first file in the provided list). 

## Convert SHAPEIT format haplotypes to VCF
     python3 consHap.py -C2V -O path_output_file_prefix -SF path_input_file_prefix
 You can use the option **-GZ** for gzipped input and output files.   
Details about SHAPEIT format are available here: https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#hapsample.  

## Convert EIGENSTRAT format to SHAPEIT format
     python3 consHap.py -C2S -O path_output_file_prefix -HF path_input_file_prefix

 You can use the option **-GZ** for gzipped input and output files.   
Details about EIGENSTRAT format are available here: https://reich.hms.harvard.edu/software/InputFileFormats.    
 

# Contact information
For any help or inquiries, please contact: ziad.albkhetan@gmail.com.
