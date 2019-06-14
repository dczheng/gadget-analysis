#!/bin/bash

ps aHux | head -1
ps aHux | grep  hybrid | grep " R"

echo "Threads Num:"  `ps aHux | grep  hybrid | grep " R" | wc -l`
