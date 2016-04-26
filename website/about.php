<!DOCTYPE html>
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
					
						
						<h3>About</h3>
						
						
						
						<div class='panel'>
							<h4>Contact Information</h4>
								<ul>
									<li>David John (john@med.uni-frankfurt.de)</li>
									<li>Dr. Shizuka Uchida (heart.lncrna@gmail.com)</li>
								</ul>
							
							<h4>Funding</h4>
								<p>This work was supported by the LOEWE Center for Cell and Gene Therapy (State of Hessen), the DFG (SFB834), and the German Center for Cardiovascular Research (DZHK).</p>
								
							<h4>Legal Statement</h4>
								<p>RnaEditor is free to use, 
								provided that the original work is properly cited. 
								It is provided "as is" without any reliability whatsoever. 
								We have taken extreme care regarding the contents that we provide in RnaEditor, 
								but if you identify a bug, please contact us. If you are commercial user, 
								please contact us: heart.lncrna@gmail.com</p>
							
							<h4>License</h4>
								<p>The MIT License (MIT)
Copyright (c) 2016 David John

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.</p>
	
				
						
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