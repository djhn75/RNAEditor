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
					
						<h1>RnaEditor</h1>
						<h3>A tool to detect editing sites from deep sequencing data</h3>
						
					<div class='panel'>
						<p>RnaEditor is an all-in-one tool to analyze RNA editing events from RNA sequencing (RNA-seq) data. 
						Raw sequencing reads (.fastq) or mapped reads (.bam) can be used as an input data. 
						RnaEditor maps the reads to the genome, calculates sequence variations, filters for “non-editing sites” and applies a cluster algorithm to detect “editing islands”. 
						</p>

					
						<div class='news'>
							<h4>News of RnaEditor</h4>
							<ul>
								<li><a href='#'>RnaEditor 1.0</a> <date>15 October 2015</date> </li>
								
							
							</ul>
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