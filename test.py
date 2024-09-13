class test1:
    def __init__(self, other_class):
        self.other_class = other_class

class test2:
    def __init__(self, some_number):
        self.some_number = some_number

class test3:
    def __init__(self):
        self.t2 = test2(5)
        self.t1 = test1(self.t2)

t3 = test3()

print(t3.t2.some_number)
print(t3.t1.other_class.some_number)

t3.t2 = None
try:
    print(t3.t2.some_number)
except:
    pass
print(t3.t1.other_class.some_number)