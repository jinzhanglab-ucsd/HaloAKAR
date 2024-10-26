// "BatchProcessFolders"
// from: https://imagej.net/ij/macros/BatchProcessFolders.txt
// This macro batch processes all the files in a folder and any
// subfolders in that folder. In this example, it runs the Subtract 
// Background command of TIFF files. For other kinds of processing,
// edit the processFile() function at the end of this macro.

   requires("1.33s"); 
   dir = getDirectory("Choose a Directory ");
   //setBatchMode(true);
   count = 0;
   countFiles(dir);
   n = 0;
   processFiles(dir);
   //print(count+" files processed");
   
   function countFiles(dir) {
      list = getFileList(dir);
      for (i=0; i<list.length; i++) {
          if (endsWith(list[i], "/"))
              countFiles(""+dir+list[i]);
          else
              count++;
      }
  }

   function processFiles(dir) {
      list = getFileList(dir);
      for (i=0; i<list.length; i++) {
          if (endsWith(list[i], "/"))
              processFiles(""+dir+list[i]);
          else {
             showProgress(n++, count);
             path = dir+list[i];
			 path2 = dir+list[i+3];
             processFile(path);
          }
      }
  }

  function processFile(path) {
       if (endsWith(path, ".001")) {
           open(path);
		   T = getTitle();
		   print(T);
		   open(path2);
		   T2 = getTitle();
		   run("Merge Channels...", "c1=["+T+"] c2=["+T2+"] create");
		   saveAs("Results", path + "_merge.tif");
           // 
           // CHANGE HERE: START
           //
        	run("ROI Manager...");
		run("Brightness/Contrast...");
       	setTool("freehand");
		
    	waitForUser("Please select ROIS. Background one first. Click OK when done");
    	run("Clear Results");
		run("Set Measurements...", "area mean min center display redirect=None decimal=3");   
		roiManager("multi-measure one");
		saveAs("Results",  path + ".csv");
		roiManager("Save", path + "_ROISet.zip");
		roiManager("Deselect");
		roiManager("Delete");
           
           // 
           // CHANGE HERE: END
           //
           close();
      }
      else{
      if (endsWith(path, ".002")) {
           open(path);
		   T = getTitle();
		   print(T);
		   open(path2);
		   T2 = getTitle();
		   run("Merge Channels...", "c1=["+T+"] c2=["+T2+"] create");
		   saveAs("Results", path + "_merge.tif");
           // 
           // CHANGE HERE: START
           //
        	run("ROI Manager...");
		run("Brightness/Contrast...");
       	setTool("freehand");
		
    	waitForUser("Please select ROIS. Click OK when done");
    	run("Clear Results");
		run("Set Measurements...", "area mean min center display redirect=None decimal=3");   
		roiManager("multi-measure one");
		saveAs("Results",  path + ".csv");
		roiManager("Save", path + "_ROISet.zip");
		roiManager("Deselect");
		roiManager("Delete");
           
           // 
           // CHANGE HERE: END
           //
           close();
      }
      else {
      if (endsWith(path, ".003")) {
           open(path);
		   T = getTitle();
		   print(T);
		   open(path2);
		   T2 = getTitle();
		   run("Merge Channels...", "c1=["+T+"] c2=["+T2+"] create");
		   saveAs("Results", path + "_merge.tif");
           // 
           // CHANGE HERE: START
           //
        	run("ROI Manager...");
		run("Brightness/Contrast...");
       	setTool("freehand");
		
    	waitForUser("Please select ROIS. Click OK when done");
    	run("Clear Results");
		run("Set Measurements...", "area mean min center display redirect=None decimal=3");   
		roiManager("multi-measure one");
		saveAs("Results",  path + ".csv");
		roiManager("Save", path + "_ROISet.zip");
		roiManager("Deselect");
		roiManager("Delete");
           
           // 
           // CHANGE HERE: END
           //
           close();
      }
      }}

  }


