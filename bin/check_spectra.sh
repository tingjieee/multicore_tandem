#!/bin/sh

echo "arg_num : $#";
echo "arg_num : $0";
echo "arg_num : $1";
echo "arg_num : $2";
cat $0 | grep PEPMASS | sort
