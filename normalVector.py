import math
import array


def factorial(r):
    if r == 0:
        return 1
    return r * factorial(r - 1)


def choose(n, k):
    # dimension = n
    # rank = k
    return factorial(n) / (factorial(k) * factorial(n - k))


def sighnPermutation(b):
    # true=even
    # false=odd
    if len(b) == 1:
        return True

    if len(b) == 2:
        if b[0] > b[1]:
            return False
        else:
            return True

    trans = 0
    for i in range(0, len(b)):
        j = i + 1
        for j in range(j, len(b)):
            if b[i] > b[j]:
                trans = trans + 1

    return (trans % 2) == 0


def findCoset(d, a):
    coset = []
    if a[0] != 1:
        for i in range(1, a[0]):
            coset.append(i)
    for i in range(len(a) - 1):
        for j in range(a[i] + 1, a[i + 1]):
            coset.append(j)
    if a[len(a) - 1] != d:
        for i in range(a[len(a) - 1] + 1, d + 1):
            coset.append(i)

    return coset


def findZero(a):
    s = a
    c = '0'
    lst = []
    for pos, char in enumerate(s):
        if (char == c):
            lst.append(pos)

    return lst


class form:
    def __init__(self, dimension, rank, coeficents):
        self.q = 0
        self.dimension = dimension
        self.rank = rank

        if len(coeficents) != choose(dimension, rank):
            raise ValueError(
                "wrong number of coeficents for a " + str(rank) + " form in " + str(
                    dimension) + " dimensions. You need " + str(
                    int(choose(self.dimension, self.rank))) + " coeficents.")
        self.coeficents = coeficents

        self.basis = [[0 for x in range(self.rank)] for y in range(int(choose(self.dimension, self.rank)))]
        one_basis = []
        for x in range(1, self.dimension + 2):
            one_basis.append(x)

        temp = [0] * (self.rank + 2)

        def combination(arr, data, start, end, index, r):
            if index == r:
                for j in range(r):
                    self.basis[self.q][j] = data[j]
                self.q += 1
                return

            i = start
            while i <= end and end - i + 1 >= r - index:
                data[index] = arr[i]
                combination(arr, data, i + 1, end, index + 1, r)
                i += 1

        combination(one_basis, temp, 0, self.dimension - 1, 0, self.rank)
        self.q = 0

    def showForm(self):
        output = ""
        for i in range(len(self.coeficents)):
            output += " " + str(self.coeficents[i]) + "("
            for j in range(self.rank):
                if j == self.rank - 1:
                    if i == len(self.coeficents) - 1:
                        output += "e_" + str(self.basis[i][j]) + ")"
                    else:
                        output += "e_" + str(self.basis[i][j]) + ") +"
                else:
                    output += "e_" + str(self.basis[i][j]) + "^"



        return output


def wedge(a, b):
    count = 0
    tensoredCoeficent = []
    tensoredBasis = []
    temp2 = [0] * (a.rank + b.rank)
    for i in range(len(a.coeficents)):
        for j in range(len(b.coeficents)):
            temp1 = a.coeficents[i] + b.coeficents[j]
            temp2 = [0] * (a.rank + b.rank)
            for k in range(len(temp2)):
                if k < a.rank:
                    temp2[k] = a.basis[i][k]

                else:
                    temp2[k] = b.basis[j][k - a.rank]

            tensoredBasis.append(temp2)
            tensoredCoeficent.append(temp1)

    dupRemovedTC = []
    dupRemovedTB = []
    repeat = True
    for i in range(len(tensoredBasis)):
        for j in range(a.rank + b.rank):
            for k in range(j + 1, a.rank + b.rank):
                if tensoredBasis[i][j] == tensoredBasis[i][k]:
                    repeat = False
        if repeat:
            dupRemovedTB.append(tensoredBasis[i])
            dupRemovedTC.append(tensoredCoeficent[i])
        repeat = True

    if a.dimension >= b.dimension:

        holder = form(a.dimension, a.rank + b.rank, [0] * int(choose(a.dimension, a.rank + b.rank)))
        outputC = [""] * int(choose(a.dimension, a.rank + b.rank))

        for i in range(int(choose(a.dimension, a.rank + b.rank))):
            for j in range(len(dupRemovedTB)):
                test = False
                count = 0
                for k in range(a.rank + b.rank):
                    for l in range(a.rank + b.rank):
                        if holder.basis[i][k] == dupRemovedTB[j][l]:
                            test = True
                    if test:
                        count += 1
                    test = False
                if count == a.rank + b.rank:
                    s = sighnPermutation(dupRemovedTB[j])
                    if s:
                        outputC[i] += "+" + dupRemovedTC[j]
                    else:
                        outputC[i] += "-" + dupRemovedTC[j]

        for i in range(len(outputC)):
            outputC[i] = "(" + outputC[i] + ")"

        return form(a.dimension, a.rank + b.rank, outputC)

    else:
        holder = form(b.dimension, a.rank + b.rank, [0] * int(choose(a.dimension, a.rank + b.rank)))
        outputC = [""] * int(choose(b.dimension, a.rank + b.rank))

        for i in range(int(choose(a.dimension, a.rank + b.rank))):
            for j in range(len(dupRemovedTB)):
                test = False
                count = 0
                for k in range(a.rank + b.rank):
                    for l in range(a.rank + b.rank):
                        if holder.basis[i][k] == dupRemovedTB[j][l]:
                            test = True
                    if test:
                        count += 1
                    test = False
                if count == a.rank + b.rank:
                    s = sighnPermutation(dupRemovedTB[j])
                    if s:
                        outputC[i] += "+" + dupRemovedTC[j]
                    else:
                        outputC[i] += "-" + dupRemovedTC[j]

            for i in range(len(outputC)):
                outputC[i] = "(" + outputC[i] + ")"

            return form(a.dimension, a.rank + b.rank, outputC)


def wedgeProduct(a):
    output = wedge(a[0], a[1])
    for i in range(2, len(a)):
        temp = output
        output = wedge(a[i], temp)
    return output


def hodge(a):
    outputC = [""] * len(a.coeficents)
    cosetBasis = [0] * len(a.coeficents)
    holder = form(a.dimension, a.dimension - a.rank, [0] * len(a.coeficents))
    for i in range(len(a.coeficents)):
        cosetBasis[i] = findCoset(a.dimension, a.basis[i])

    for i in range(len(a.coeficents)):
        for j in range(len(a.coeficents)):
            if cosetBasis[i] == holder.basis[j]:
                temp = a.basis[i]
                temp2 = temp + cosetBasis[i]
                if sighnPermutation(temp2):
                    outputC[j] = "+(" + a.coeficents[i] + ")"
                if not (sighnPermutation(temp2)):
                    outputC[j] = "-(" + a.coeficents[i] + ")"
                temp = []
    return form(a.dimension, a.dimension - a.rank, outputC)


def nDTest(n):
    coeficents = [[0] * n] * (n - 2)
    forms = [0] * (n - 2)
    for i in range(1, n - 1):
        for j in range(1, n):
            coeficents[i - 1][j - 1] = "(c_" + str(j) + ";" + str(i) + ")"

    for i in range(n - 2):
        coeficents[i][n - 1] = "0"

    for i in range(n - 2):
        forms[i] = form(n, 1, coeficents[i])

    sformc = [0] * (n - 1)
    for i in range(1, n):
        sformc[i - 1] = "(c_" + str(i) + "-p_" + str(i) + ")"

    sformc.append("(-p_" + str(n) + ")")
    sform = form(n, 1, sformc)

    F = []

    for i in range(n):
        F.append("(x_"+str(i+1)+")")

    return wedge(wedge(wedgeProduct(forms), sform),form(n,1,F))


print(nDTest(488).showForm())

