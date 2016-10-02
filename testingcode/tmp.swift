import files;

files = glob("/home/azza/swift-project/Swift-Variant-Calling/*.swift");

size(files);

trace ("#########\t\t" + size(files) + "#########");


trace ("#########\t\t" , filename(files[0]) , "#########\n\n");
