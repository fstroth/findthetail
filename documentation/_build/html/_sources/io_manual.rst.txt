Input/Output
============

Input
-----

The input has to be a 1-dimensional time series. The data has to in the form of a numpy array. There is no need for any preperation of the data in the form of sorting or making sure no value is given more than twice. All this prepatrions will be done by the FTT procedure.

The preperation procedure of consits of two steps:

    1. Check if any value exists more than once, if so add small random number below the significant digit of the data
    2. Sort the data in descending order

Output
------

From the information of the analysis a html report is generated which contains all important information. For a more detailed infomation on the report see the next section.