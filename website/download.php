<html>
	<?php
        /* Write Header which includes all scripts and Stylesheets */
        include ('html/header.html');
    ?>
	
	<body>
	<img src="img/rnaEditor_64x64.png" style="position: absolute;" id="waitImg" width=50"/>
       <div class="row">
       	<div class="small-12 large-12 columns">
            <?php
                /* Write Menu */
                include ('html/menu.html');
            ?>
		 </div>
		</div>
		
		<div class="row">
        	<div class="small-12 large-12 columns">    
            
	            <div class="content">
					
						
						<h3>Download indexes and annotations</h3>
						
						<div class='panel'>
							
							<p>RnaEditor requires a set of annotation files and databases to detect editing sites.
							Either download one of the following annotation bundles or download the files manually by executing the commands below.</p>
							
							<h5> Human </h5>
								<ul>
									<li><a href="http://141.2.194.197/rnaeditor_annotations/GRCH38.tar.gz">GRCH38</a></li>
									<li><a href="http://141.2.194.197/rnaeditor_annotations/GRCH37.tar.gz">GRCH37</a></li>
								</ul>
							<h5> Mouse </h5>
								<ul>
									<li><a href="http://141.2.194.197/rnaeditor_annotations/GRCM38.tar.gz">GRCHM38</a></li>
									
								</ul>
							<br>
							<h5>Unix commands to download GRCH38 manually</h5>
								<li>Reference Genome:</li>
									<div class="codebox">
									    <code>   
									        wget -qO- ftp://ftp.ensembl.org/pub/release-83/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz |gunzip -c > Homo_sapiens.GRCh38.dna.primary_assembly.fa
									    bwa index Homo_sapiens.GRCh38.dna.primary_assembly.fa; samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
									    
									    </code>
									</div>
								<li>Gene Annotation:</li>
									<div class="codebox">
									    <code>   
									       wget -qO- ftp://ftp.ensembl.org/pub/release-83/gtf/homo_sapiens/Homo_sapiens.GRCh38.83.gtf.gz |gunzip -c > Homo_sapiens.GRCh38.83.gtf
										</code>
									</div>
								<li>dbSNP Database:</li>
									<div class="codebox">
									    <code>   
									    	wget -qO- ftp://ftp.ensembl.org/pub/release-83/variation/vcf/homo_sapiens/Homo_sapiens.vcf.gz |gunzip -c |awk 'BEGIN{FS="\t";OFS="\t"};match($5,/\./){gsub(/\./,"N",$5)};$5 == "" && $1 !~ /^#/ {gsub("","N",$5)};$3 ~ /rs193922900/ {$5="TN"};$3 ~ /rs59736472/ {$5="AN"};$5 ~ /H/ {gsub(/H/,"N",$5)};{print $0}' dbSNP.vcf   
										</code>
									</div>
								<li>ESP Database:</li>
									<div class="codebox">
									    <code>   
									       wget -qO- ftp://ftp.ensembl.org/pub/release-83/variation/vcf/homo_sapiens/ESP65*.vcf.gz |gunzip -c |grep -v ^## |grep -v rs[0-9][0-9] > ESP.vcf
										</code>
									</div>
								<li>HAPMAP Database:</li>
									<div class="codebox">
									    <code>   
									       wget -qO- ftp://ftp.ensembl.org/pub/release-83/variation/vcf/homo_sapiens/*HAPMAP*.vcf.gz |gunzip -c |grep -v ^## |grep -v rs[0-9][0-9] > HAPMAP.vcf
										</code>
									</div>
								<li>Download Alu Regions:</li>
									
									    <ul>
									    <li>#download Alu regions from the repeat masker
											<ul><li>Link: http://genome.ucsc.edu/cgi-bin/hgTables</li>
												<li>group: Variation and Repeats</li>
												<li>track: RepeatMasker</li>
												<li>table: rmsk</li>
												<li>output format: BED</li>
											</ul>
									    
									    </ul>
									<div class="codebox">
										<code> 
											#run this awk command to make the alu region compatible to the ucsc annotation
											<br>
											awk 'BEGIN{FS="\t";OFS="\t"} match($1,/chr/){$1 = substr($1,4)}{print $0}' yourFile.bed
										</code>
									</div>
						
						</div>
	            </div>
			</div>
		</div>


		<?php
                /* Write Menu */
                include ('html/footer.html');
            ?>
	</body>


</html>