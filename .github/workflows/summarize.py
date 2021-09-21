import sys
input_log = sys.argv[1]
output = sys.argv[2]

log = open(input_log).read().split("\n")
summary = {}
for l in log:
	if l.startswith("*"): continue
	items = l.split(":")
	if len(items) < 5: continue
	if not items[3] in summary: summary[items[3]] = {"text": items[4], "count": 0}
	summary[items[3]]["count"] += 1

import json
json.dump(summary, open(output, "w"), sort_keys=True, indent=4, separators=(',', ': '))
