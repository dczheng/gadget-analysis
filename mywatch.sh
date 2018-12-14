#!/usr/bin/fish

ps aHux | head -1
ps aHux | grep  gadget-analysis | grep " R"
set s (ps aHux | grep  gadget-analysis | grep " R" | wc -l)

echo "Threads Num:" $s
