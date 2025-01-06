// Batch count nuclei from a directory

// Define functions ---------

// Define function for cropping 928 square pixel images
function crop_image(raw_image_dir, cropped_image_dir, filename) {
        open(raw_image_dir + filename);
        makeRectangle(50, 50, 928, 928);
        run("Crop");
        saveAs("tif", cropped_image_dir + filename);
        close();
}

// Define function for processing and counting nuclei
function process_and_count_nuclei(cropped_image_dir, processed_image_dir, filename) { 
	open(cropped_image_dir + filename);
	setOption("BlackBackground", true);
	run("Subtract Background...", "rolling=12");
	setAutoThreshold("Default dark no-reset");
	//setThreshold(20, 255);
	run("Convert to Mask");
	run("Fill Holes");
	run("Convert to Mask");
	run("Watershed");
	run("Analyze Particles...", "size=3-50 circularity=0.50-1.00 show=Overlay display exclude clear summarize");
	saveAs("tif", processed_image_dir + filename);
    close();
}

// Define image directories to be used -----
raw_image_dir = "/Users/tsoleary/OneDrive - University of Vermont/Research/HEATER/nuclei_concentration_10X_samples/raw_images/";
cropped_image_dir = "/Users/tsoleary/OneDrive - University of Vermont/Research/HEATER/nuclei_concentration_10X_samples/cropped_images/";
processed_image_dir = "/Users/tsoleary/OneDrive - University of Vermont/Research/HEATER/nuclei_concentration_10X_samples/processed_images/";

// Crop all images in raw_image_dir -----
list = getFileList(raw_image_dir);

for (i = 0; i < list.length; i++){
	crop_image(raw_image_dir, cropped_image_dir, list[i]);
}

// Process all images in cropped image dir ----
list = getFileList(cropped_image_dir);
for (i = 0; i < list. length; i++){
	process_and_count_nuclei(cropped_image_dir, processed_image_dir, list[i]);
}
