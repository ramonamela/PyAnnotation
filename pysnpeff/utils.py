#!/usr/bin/env python3


def _annotate_streams(input_stream, output_stream, config):
    """Annotates the input_stream to the output_stream following the config
    directives.
    Parameters
    ----------
    input_stream : TextIOBase
        Input stream.
    output_stream : TextIOBase
        Output stream.
    config : dict
        Dictionary describing the annotation step to be performed
    """
    pass

def annotate(input_file, output_file, json_file):
    """Annotates the input_file to the output_file following the configuration
    present in json_file.
    Parameters
    ----------
    input_stream : TextIOBase
        Input stream.
    output_stream : TextIOBase
        Input stream.
    config : dict
        Dictionary containing the necessary information to perform the
        annotation
    """
    pass

def main(*args, **kwargs):
    pass

if __name__ == "__main__":
    ## Check input parameters and launch the main function
    #main()