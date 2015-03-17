# SPADE Analysis with d3.js Visualization
THe following project outlines how to run a quick SPADE analysis using the [flowCore](http://www.bioconductor.org/packages/release/bioc/html/flowCore.html) and [SPADE](http://www.bioconductor.org/packages/release/bioc/html/spade.html) pacakges

# Procedure
 1. download files to analyze
 2. run SPADE analysis
 3. summarize SPADE information and save as JSON
 4. visualize

# Run the web app
Using python3 on windows, from the console go to the /web folder and run the following command

```
python3 -m http.server
```

then point your browser to 

```
http://localhost:8000/spade.html
```

# Acknolwedgments
Mike Bostock for the [d3.js package](https://github.com/mbostock/d3/wiki/Gallery)
Bernd Bodenmiller for his [fantastic publication](http://www.nature.com/nbt/journal/v30/n9/full/nbt.2317.html)