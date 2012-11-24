temp = """
def get_%s(data) :
    return data[%s:%s+%s]
"""

def foo(name, offset, leng):
    exec temp % (name, offset, offset, leng) in globals()

t = foo('a',4,4)
print t
