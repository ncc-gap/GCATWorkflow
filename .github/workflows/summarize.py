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
print("""The message type can be:
[I]nformational messages that Pylint emits (do not contribute to your analysis score)
[R]efactor for a "good practice" metric violation
[C]onvention for coding standard violation
[W]arning for stylistic problems, or minor programming issues
[E]rror for important programming issues (i.e. most probably bug)
[F]atal for errors which prevented further processing""")
