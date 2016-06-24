RNAEditor is tool to analyze RNA editing events from RNA-seq data.
It accepts **Fastq** and **Bam** files as input, handles all the analysis and filter steps for you and outputs comprehensive statistics and **editing Islands**.


#### 1.) Reqiurements

Install following tools to **/usr/local/bin/**:

+ BWA: [Download](https://sourceforge.net/projects/bio-bwa/files/latest/download)
+ Picard Tools: [Donwload](https://sourceforge.net/projects/picard/files/picard-tools/1.119/picard-tools-1.119.zip/download)
	**Move all .jar files to /usr/local/bin/picard-tools/** (use version <= 1.119)
+ GATK: [Donwload](https://www.broadinstitute.org/gatk/download/auth?package=GATK)
+ Blat:
	+ Binarys:
		+ [MacOSX_x86_64](http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/blat/blat)
		+ [Linux.x86_64](http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/blat)
	+ Source:
		+ [blatSrc35.zip](https://users.soe.ucsc.edu/~kent/src/blatSrc35.zip)
+ Bedtools: [Download]( http://bedtools.readthedocs.io/en/latest/content/installation.html)

Also you need to make sure that you installed **pysam, python-qt4, matplotlib, numpy**.
If pip is available just run the following comannds to fulfill the requirements.
`$ pip instal pysam`
`$ sudo apt-get install python-qt4`
`$ sudo  matplotlib` (if this fails run `sudo apt-get install libfreetype6-dev libxft-dev` first)
`$ pip install pysam`


#### 2.)Install
To install RNAEditor either donwload the prebuild Application bundle or clone git repository
+	MacOSX: [RNAEditor.app](http://rnaeditor.uni-frankfurt.de/src/RNAEditor.dmg)
+	Linux: [RNAEditor.tar.gz](http://rnaeditor.uni-frankfurt.de/src/RNAEditor.tar.gz)
+	Source: [GitHub]: (https://github.com/djhn75/RNAEditor)

### 3.) Start RNAEditor

+ run`python RNAEditor.py`

#### 4.) Run Analysis
To detect editing from you NGS data simply drop all your FASTQ-Files to the RNAEditor window set the parameters and press START


#### 5.) Legal Statement
RNAEditor is free to use without login information, provided that the original work is properly cited.
It is provided “as is” without any reliability whatsoever.
We have taken extreme care regarding the contents that we provide in RNAEditor, but if you identify a bug, please contact us.
If you are commercial user, please contact us: uchida@med.uni-frankfurt.de

### 6.) LICENSE

The MIT License

Copyright (c) 2016 David John

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.