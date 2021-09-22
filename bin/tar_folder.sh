#!/bin/bash

destination=$1
folder=$2

echo -e "[$(date)\ttaring $folder to $destination"
tar -zcf $destination $folder

if [ "$?" -ne 0 ]
then
	echo "Command returned non-zero return value. File should be checked. $destination"
fi

