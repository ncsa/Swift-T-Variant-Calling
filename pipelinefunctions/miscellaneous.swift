// These are functions necessary for the pipeline, but hard to classify into a job.

// Reading the runfile parameters:
(string data[string]) getConfigVariables(string lines[])
{
        foreach line in lines
        {
                string keyValuePair[] = split(line, "=");
                string name = keyValuePair[0];
                string value = keyValuePair[1];
                data[name] = value;
        }
}


app (file output) find_files (string dir, string pattern){
        "/usr/bin/find" dir "-name" pattern @stdout=output;
}

