<html>
	<?php
        /* Write Header which includes all scripts and Stylesheets */
        include ('html/header.html');
    ?>
    <style>

      img, p, u { margin-left: 25px; }
      h4 {margin-top: 60px;}
    </style>
	
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
					
						
						<h3>RNAEditor Documentation</h3>
						
						<div class='panel'>
							
								
						  <h3>Output Files</h3>
                                 <h4> VCF File:</h4>
                                    <p>A standard vcf file with all the editing sites for further analysis.</p>
                                    <u>Example: Sample.vcf</u>
                                        <img src="img/vcf-file.png"/>
                                <h4> GVF File:</h4>
                                    <p><b>G</b>ene <b>V</b>ariation <b>F</b>ile hold additional informations for each editing site, like  gene names, segments, number of total reads, number of edited reads and the ratio of editing.</p>
                                    <u>Example: Sample.gvf</u>
                                        <img src="img/gvf-file.png"/>
                                <h4> Summary File:</h4>
                                    <p>Shows the sum of editing sites per segment for each gene.</p>
                                    <u>Example: Sample.summary</u><br>
                                        <img src="img/summary-file.png"/>
                                <h4> Editing Island File:</h4>
                                    <p></p>
                                    <u>Example: Sample.clusters.bed</u>
						              <img src="img/island-file.png"/>
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