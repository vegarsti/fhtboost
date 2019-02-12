import sys

print("c(" + ", ".join(line.strip() for line in sys.stdin) + ")")
