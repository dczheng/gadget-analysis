#!/bin/bash

ps aHux | head -1
ps aHux | grep  gadget-analysis | grep " R"

echo "Threads Num:"  `ps aHux | grep  gadget-analysis | grep " R" | wc -l`
