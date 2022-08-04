import pickle

class A(object):
    def __init__(self, important_data):
        self.important_data = important_data

        # Add data which cannot be pickled:
        self.func = lambda: 7

        # Add data which should never be pickled, because it expires quickly:
        self.is_up_to_date = False

    def add_two(self, x):
        return x+2

    def __getstate__(self):
        return [self.important_data] # only this is needed

    def __setstate__(self, state):
        self.important_data = state[0]

        self.func = lambda: 7  # just some hard-coded unpicklable function
        self.func_2 = self.add_two

        self.is_up_to_date = False  # even if it was before pickling


if __name__ == "__main__":
    a1 = A('very important')
    with open('s.dat', 'wb') as outf:
        pickle.dump(a1, outf)

    # s = pickle.dumps(a1)  # calls a1.__getstate__()
    # a2 = pickle.loads(s)  # calls a1.__setstate__(['very important'])

    infile = open('s.dat', 'rb')
    a2 = pickle.load(infile)
    infile.close()
    print(a2)
    print(a2.important_data)
    print(a2.func())
    print(a2.add_two(10))
