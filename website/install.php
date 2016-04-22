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
					
						
						<h3>Installation</h3>
						
						<div class='panel'>
							<p>To install RnaEditor follow these steps :</p>
							
							<h4>Required Software:</h4>
								<p>RnaEditor will ask for a source directly where the following software products have to be installed:</p>
								<ul>
									<li><a href="http://bio-bwa.sourceforge.net">BWA</a></li>
									<li><a href="http://broadinstitute.github.io/picard/">Picard Tools</a></li>
									<li><a href="https://www.broadinstitute.org/gatk/">GATK</a></li>
									<li><a href="http://hgdownload.soe.ucsc.edu/downloads.html#utilities_downloads">Blat</a></li>
									<li><a href="http://bedtools.readthedocs.org/en/latest/">Bedtools</a></li>										
								</ul>

									
								
							<h4>RnaEditor:</h4>
								<ul>
									<li>Unix:</li>
									<li>Mac:</li>
									<li>Windows:</li>
								
								
								</ul>
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