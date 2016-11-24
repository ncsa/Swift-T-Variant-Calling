import sys

print("")
print("ARGUMENT CHECKER")
print("argv: " + str(sys.argv))
n = len(sys.argv)
print("arg count: " + str(n))
for i in range(1,n):
    a = sys.argv[i]
    print "argument: " + str(i)      + \
         " length: "   + str(len(a)) + \
         " value: '"   + a + "'"
