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
            
	    <!-- Table of Content-->
        <div id='toc' class='toc'>
        <div id='toctitle'><h2>Content</h2></div>
        <ul>
            <li class='toclevel-1 tocsection-1'><a href='#guiMode'><span class='tocnumber'>1</span> <span class='toctext'>Interactive mode</span></a></li>
            <li class='toclevel-1 tocsection-2'><a href='#consoleMode'><span class='tocnumber'>2</span> <span class='toctext'>Console mode</span></a></li>
                        
                    </div>	            
                    
                    
                    
                    
                    
                    <div class="content">
					
						<p> RnaEditor can be used in two ways. Either you start it with a user inteface or you use it from the command line.
						Both methods will be explained below.</p>
						
						
						<h3 id='guiMode'>Interactive mode</h3>
							<div class='panel'>
								<p>To start RnaEitor with a user interface just run <code>python RnaEditor.py</code> or double click the program symbol on Mac or Windows.</p>
								
								<p>Follow these steps to run you analisys:</p>
								<ul>
									<li>1 Download your <a href="download.php"><b>annotation bundle</b></a> and extract it</li>
									<li>2 Set the correct paths in the file <b>configuration.txt</b> inside of your bundle</li>
									<li>3 Drop <b>configuration.txt</b> inside of RnaEditor</li>
									<li>4 Drop your Fastq files inside of RnaEditor</li>
									<li>5 Check the parameters again and Press start then wait for the results</li>
								<ul>
						
							</div>
						
						
						<h3 id='consoleMode'>Console mode</h3>
							<div class='panel'>
								<p>For an overview on how to run RnaEditor just run <code>python RnaEditor.py -h</code>. This will print the following help page. </p>								
								<div class="codebox">

								<code>
								usage: RnaEditor [-h] -i Fastq-Files [Fastq-Files ...] -c Configuration File<br><br>

									 RnaEditor: easily detect editing sites from deep sequencing data'<br>
									----------------------------------------------------------------<br>
									    run without arguments to start the user interface.<br>
									<br>
									optional arguments:<br><br>
									  -h, --help            show this help message and exit<br>
									  -i Fastq-Files [Fastq-Files ...], --input Fastq-Files [Fastq-Files ...]
									                        Input fastq files (maximum two for paire-end-
									                        sequencing)<br>
									  -c Configuration-File, --conf Configuration-File
									                        Configuration File used to read Parameters for
									                        RnaEditor<br>							
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