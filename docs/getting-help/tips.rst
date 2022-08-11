.. _getting_help-faqs::

===============
Tips and Tricks
===============

On this page you can find short answers to frequently asked questions about phylociraptor.

----------------------------------------------
Controlling the number of threads in -t serial
----------------------------------------------

By default phylociraptor will use the maximum number of available threads when it is run in :bash:`-t serial` mode. If you would like to change this you can use :bash:`-t serial=4` to restrict it to use only four threads for eaxmple. This only applies to parts which alread use multithreading, such as the tree calculation steps.

