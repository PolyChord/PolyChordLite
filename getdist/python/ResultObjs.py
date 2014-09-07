import paramNames, decimal

class textFile:

    def __init__(self, lines=[]):
        self.lines = lines

    def write(self, outfile):
        textFileHandle = open(outfile, 'w')
        textFileHandle.write("\n".join(self.lines))
        textFileHandle.close()

def texEscapeText(string):
        return string.replace('_', '{\\textunderscore}')


def float_to_decimal(f):
    # http://docs.python.org/library/decimal.html#decimal-faq
    "Convert a floating point number to a Decimal with no loss of information"
    n, d = f.as_integer_ratio()
    numerator, denominator = decimal.Decimal(n), decimal.Decimal(d)
    ctx = decimal.Context(prec=60)
    result = ctx.divide(numerator, denominator)
    while ctx.flags[decimal.Inexact]:
        ctx.flags[decimal.Inexact] = False
        ctx.prec *= 2
        result = ctx.divide(numerator, denominator)
    return result

def numberFigs(number, sigfig):
    # http://stackoverflow.com/questions/2663612/nicely-representing-a-floating-point-number-in-python/2663623#2663623
    assert(sigfig > 0)
    try:
        d = decimal.Decimal(number)
    except TypeError:
        d = float_to_decimal(float(number))
    sign, digits = d.as_tuple()[0:2]
    if len(digits) < sigfig:
        digits = list(digits)
        digits.extend([0] * (sigfig - len(digits)))
    shift = d.adjusted()
    result = int(''.join(map(str, digits[:sigfig])))
    # Round the result
    if len(digits) > sigfig and digits[sigfig] >= 5: result += 1
    result = list(str(result))
    # Rounding can change the length of result
    # If so, adjust shift
    shift += len(result) - sigfig
    # reset len of result to sigfig
    result = result[:sigfig]
    if shift >= sigfig - 1:
        # Tack more zeros on the end
        result += ['0'] * (shift - sigfig + 1)
    elif 0 <= shift:
        # Place the decimal point in between digits
        result.insert(shift + 1, '.')
    else:
        # Tack zeros on the front
        assert(shift < 0)
        result = ['0.'] + ['0'] * (-shift - 1) + result
    if sign:
        result.insert(0, '-')
    return ''.join(result)

class numberFormatter():
    def __init__(self, sig_figs=4):
        self.sig_figs = sig_figs
        self.separate_limit_tol = 0.1

    def namesigFigs(self, value, limplus, limminus, wantSign=True):
        frac = limplus / (abs(value) + limplus)
        err_sf = 2
        if value >= 20 and frac > 0.1 and limplus >= 2: err_sf = 1

        plus_str = self.formatNumber(limplus, err_sf, wantSign)
        minus_str = self.formatNumber(limminus, err_sf, wantSign)
        sf = self.sig_figs
        if frac > 0.1 and value < 100 and value >= 20: sf = 2
        elif frac > 0.01 and value < 1000: sf = 3
#        if abs(value) < 1 and limplus - limminus > abs(value): sf = 2
        res = self.formatNumber(value, sf)
        maxdp = max(self.decimal_places(plus_str), self.decimal_places(minus_str))
        while abs(value) < 1 and maxdp < self.decimal_places(res):
            sf -= 1
            if sf == 0:
                res = ('%.' + str(maxdp) + 'f') % value
                if (float(res) == 0.0): res = ('%.' + str(maxdp) + 'f') % 0
                break
            else: res = self.formatNumber(value, sf)

        while self.decimal_places(plus_str) > self.decimal_places(res):
            sf += 1
            res = self.formatNumber(value, sf)
        return (res, plus_str, minus_str)

    def formatNumber(self, value, sig_figs=None, wantSign=False):
        if sig_figs is None:
            sf = self.sig_figs
        else: sf = sig_figs
        s = numberFigs(value, sf)
        if wantSign:
            if s[0] != '-' and float(s) < 0: s = '-' + s
            if  float(s) > 0: s = '+' + s
        return s

    def decimal_places(self, s):
        i = s.find('.')
        if i > 0: return len(s) - i - 1
        return 0

    def plusMinusLimit(self, limit, upper, lower):
        return limit != 1 or abs(abs(upper / lower) - 1) > self.separate_limit_tol


class tableFormatter(object):
    def __init__(self):
        self.border = '|'
        self.endofrow = '\\\\'
        self.hline = '\\hline'
        self.paramText = 'Parameter'
        self.aboveTitles = self.hline
        self.majorDividor = '|'
        self.minorDividor = '|'
        self.colDividor = '||'
        self.belowTitles = ''
        self.headerWrapper = " %s"
        self.noConstraint = '---'
        self.spacer = ' '  # just to make output more readable
        self.colSeparator = self.spacer + '&' + self.spacer
        self.numberFormatter = numberFormatter()

    def getLine(self, position=None):
        if position is not None and hasattr(self, position): return getattr(self, position)
        return self.hline

    def belowTitleLine(self, colsPerParam, numResults=None):
        return self.getLine("belowTitles")

    def startTable(self, ncol, colsPerResult, numResults):
        part = self.majorDividor + (" c" + self.minorDividor) * (colsPerResult - 1) + ' c'
        return '\\begin{tabular} {' + self.border + " l " + part * numResults + (self.colDividor + " l " + part * numResults) * (ncol - 1) + self.border + '}'

    def endTable(self):
        return '\\end{tabular}'

    def titleSubColumn(self, colsPerResult, title):
        return ' \\multicolumn{' + str(colsPerResult) + '}{' + self.majorDividor + 'c' + self.majorDividor + '}{' + self.formatTitle(title) + '}'

    def formatTitle(self, title):
        return '\\bf ' + texEscapeText(title)

    def textAsColumn(self, txt, latex=False, separator=False, bold=False):
        wid = len(txt)
        if latex:
            wid += 2
            if bold: wid += 11
        res = txt + self.spacer * max(0, 28 - wid)
        if latex:
            if bold: res = '{\\boldmath$' + res + '$}'
            else:  res = '$' + res + '$'
        if separator: res += self.colSeparator
        return res

class openTableFormatter(tableFormatter):

    def __init__(self):
        tableFormatter.__init__(self)
        self.border = ''
        self.aboveTitles = r'\noalign{\vskip 3pt}' + self.hline + r'\noalign{\vskip 1.5pt}' + self.hline + r'\noalign{\vskip 5pt}'
        self.belowTitles = r'\noalign{\vskip 3pt}' + self.hline
        self.aboveHeader = ''
        self.belowHeader = self.hline
        self.minorDividor = ''
        self.belowFinalRow = ''

    def titleSubColumn(self, colsPerResult, title):
        return ' \\multicolumn{' + str(colsPerResult) + '}{' + 'c' + '}{' + self.formatTitle(title) + '}'

class noLineTableFormatter(openTableFormatter):

    def __init__(self):
        openTableFormatter.__init__(self)
        self.aboveHeader = ''
#        self.belowHeader = r'\noalign{\vskip 5pt}'
        self.minorDividor = ''
        self.majorDividor = ''
        self.belowFinalRow = self.hline
        self.belowBlockRow = self.hline
        self.colDividor = '|'
        self.hline = ''

    def belowTitleLine(self, colsPerParam, numResults=None):
        return r'\noalign{\vskip 3pt}\cline{2-' + str(colsPerParam * numResults + 1) + r'}\noalign{\vskip 3pt}'


class resultTable():

    def __init__(self, ncol, results, limit=2, tableParamNames=None, titles=None, formatter=None,
                 numFormatter=None, blockEndParams=None, paramList=None, refResults=None):
# results is a margeStats or bestFit table
        self.lines = []
        if formatter is None: self.format = noLineTableFormatter()
        else: self.format = formatter
        self.ncol = ncol
        if tableParamNames is None:
            self.tableParamNames = results[0]
        else: self.tableParamNames = tableParamNames
        if paramList is not None: self.tableParamNames = self.tableParamNames.filteredCopy(paramList)
        if numFormatter is not None: self.format.numFormatter = numFormatter

        self.results = results
        self.boldBaseParameters = True
        self.colsPerResult = len(results[0].getColumnLabels(limit))
        self.colsPerParam = len(results) * self.colsPerResult
        self.limit = limit
        self.refResults = refResults

        nparams = self.tableParamNames.numParams()
        numrow = nparams / ncol
        if nparams % ncol != 0: numrow += 1
        rows = []
        for par in self.tableParamNames.names[0:numrow]:
            rows.append([par])
        for col in range(1, ncol):
            for i in range(numrow * col, min(numrow * (col + 1), nparams)):
                rows[i - numrow * col].append(self.tableParamNames.names[i]);

        self.lines.append(self.format.startTable(ncol, self.colsPerResult, len(results)))
        if titles is not None: self.addTitlesRow(titles)
        self.addHeaderRow()
        for row in rows[:-1]:
            self.addFullTableRow(row)
            if ncol == 1 and blockEndParams is not None and row[0].name in blockEndParams: self.addLine("belowBlockRow")
            else: self.addLine("belowRow")
        self.addFullTableRow(rows[-1])
        self.addLine("belowFinalRow")
        self.endTable()


    def addFullTableRow(self, row):
        txt = self.format.colSeparator.join(self.paramLabelColumn(param) + self.paramResultsTex(param) for param in row)
        if not self.ncol == len(row):
            txt += self.format.colSeparator * ((1 + self.colsPerParam) * (self.ncol - len(row)))
        self.lines.append(txt + self.format.endofrow)

    def addLine(self, position):
        if self.format.getLine(position) is None:  # no line is appended if the attribute is None
            return self.lines
        else:
            return self.lines.append(self.format.getLine(position))

    def addTitlesRow(self, titles):
        self.addLine("aboveTitles")
        cols = [self.format.titleSubColumn(1, '')]
        cols += [self.format.titleSubColumn(self.colsPerResult, title) for title in titles]
        self.lines.append(self.format.colSeparator.join(cols * self.ncol) + self.format.endofrow)

#        res = self.format.titleSubColumn(1, '') + self.format.colSeparator + self.format.colSeparator.join(self.format.titleSubColumn(self.colsPerResult, title) for title in titles)
#        self.lines.append(((self.format.colSeparator + res) * self.ncol)[1:] + self.format.endofrow)
        belowTitleLine = self.format.belowTitleLine(self.colsPerResult, self.colsPerParam / self.colsPerResult)
        if belowTitleLine:
            self.lines.append(belowTitleLine)

    def addHeaderRow(self):
        self.addLine("aboveHeader")
        cols = [self.format.headerWrapper % self.format.paramText]
        for result in self.results:
            cols += [self.format.headerWrapper % s for s in result.getColumnLabels(self.limit)]
        self.lines.append(self.format.colSeparator.join(cols * self.ncol) + self.format.endofrow)

#        res = self.format.colSeparator + self.format.headerWrapper % self.format.paramText
#        for result in self.results:
#            res += self.format.colSeparator + self.format.colSeparator.join([self.format.headerWrapper % s for s in result.getColumnLabels(self.limit)])
#        self.lines.append((res * self.ncol).replace(self.format.colSeparator,'',1) + self.format.endofrow)
        self.addLine("belowHeader")

    def paramResultsTex(self, param):
        return self.format.colSeparator.join(self.paramResultTex(result, param) for result in self.results)

    def paramResultTex(self, result, p):
        values = result.texValues(self.format, p, self.limit, self.refResults)
        if values is not None:
            if len(values) > 1: txt = self.format.textAsColumn(values[1], True, separator=True)
            else: txt = ''
            txt += self.format.textAsColumn(values[0], values[0] != self.format.noConstraint)
            return txt
        else: return self.format.textAsColumn('') * len(result.getColumnLabels(self.limit))

    def paramLabelColumn(self, param):
        return  self.format.textAsColumn(param.label, True, separator=True, bold=not param.isDerived)

    def endTable(self):
        self.lines.append(self.format.endTable())

    def tableTex(self):
        return "\n".join(self.lines)

    def writeTable(self, fname):
        textFile(self.lines).write(fname)


class paramResults(paramNames.paramList): pass

class bestFit(paramResults):

    def __init__(self, fileName=None, setParamNameFile=None, want_fixed=False):
        paramResults.__init__(self)
        if fileName is not None: self.loadFromFile(fileName, want_fixed=want_fixed)
        if setParamNameFile is not None: self.setLabelsAndDerivedFromParamNames(setParamNameFile)

    def getColumnLabels(self, limit=None):
        return ['Best fit']

    def loadFromFile(self, filename, want_fixed=False):
        textFileLines = self.fileList(filename)
        first = textFileLines[0].strip().split('=')
        if first[0].strip() == 'weight':
            self.weight = float(first[1].strip())
            del(textFileLines[0])
            first = textFileLines[0].strip().split('=')
        if first[0].strip() != '-log(Like)': raise Exception('Error in format of parameter (best fit) file')
        self.logLike = float(first[1].strip())
        isFixed = False
        isDerived = False
        self.chiSquareds = []
        chunks = 0
        if len(textFileLines[1].strip()) > 0: del(textFileLines[1])  # if it has chi2 line as well
        for ix in range(2, len(textFileLines)):
            line = textFileLines[ix]
            if len(line.strip()) == 0:
                chunks += 1
                isFixed = not isFixed
                isDerived = True
                if chunks == 3:
                    if ix + 2 >= len(textFileLines): break
                    for likePart in textFileLines[ix + 2:]:
                        if len(likePart.strip()) != 0:
                            (chisq, name) = [s.strip() for s in likePart.split(None, 2)][1:]
                            name = [s.strip() for s in name.split(':', 1)]
                            if len(name) > 1:
                                (kind, name) = name
                            else: kind = ''
                            name = name.replace('.clik', '').replace('.cldf', '')
                            self.chiSquareds.append((kind, name, float(chisq)))
                    break
                continue
            if not isFixed or want_fixed:
                param = paramNames.paramInfo()
                param.isDerived = isDerived
                (param.number, param.best_fit, param.name, param.label) = [s.strip() for s in line.split(None, 3)]
                param.number = int(param.number)
                param.best_fit = float(param.best_fit)
                self.names.append(param)

    def sortedChiSquareds(self):
        likes = dict()
        for (kind, name, chisq) in self.chiSquareds:
            if not kind in likes: likes[kind] = []
            likes[kind].append((name, chisq))
        return sorted(likes.iteritems())

    def chiSquareForKindName(self, kind, name):
        for (akind, aname, chisq) in self.chiSquareds:
            if akind == kind and aname == name: return chisq
        return None


    def texValues(self, formatter, p, limit=None, refResults=None):
        param = self.parWithName(p.name)
        if param is not None: return [formatter.numberFormatter.formatNumber(param.best_fit)]
        else: return None

class paramLimit():
    def __init__(self, minmax, tag='two'):
        self.lower = minmax[0]
        self.upper = minmax[1]
        self.twotail = tag == 'two'
        self.onetail_upper = tag == '>'
        self.onetail_lower = tag == '<'

class margeStats(paramResults):

    def loadFromFile(self, filename):
        print filename
        textFileLines = self.fileList(filename)
        lims = textFileLines[0].split(':')[1]
        self.limits = [float(s.strip()) for s in lims.split(';')]
        self.hasBestFit = False
        for line in textFileLines[3:]:
            if len(line.strip()) == 0: break
            param = paramNames.paramInfo()
            items = [s.strip() for s in line.split(None, len(self.limits) * 3 + 3)]
            param.name = items[0]
            param.mean = float(items[1])
            param.err = float(items[2])
            param.label = items[-1]
            param.limits = []
            for i in range(len(self.limits)):
                param.limits.append(paramLimit([float(s) for s in items[3 + i * 3:5 + i * 3] ], items[5 + i * 3]))
            self.names.append(param)


    def addBestFit(self, bf):
        self.hasBestFit = True
        self.logLike = bf.logLike
# the next line deletes parameters not in best-fit; this is good e.g. to get rid of yhe from importance sampled result
        self.names = [x for x in self.names if bf.parWithName(x.name) is not None]
        for par in self.names:
            param = bf.parWithName(par.name)
            par.best_fit = param.best_fit
            par.isDerived = param.isDerived

    def getColumnLabels(self, limit=2):
        if self.hasBestFit: res = ['Best fit']
        else: res = []
        number_string = str(round(float(self.limits[limit - 1]) * 100))
        if number_string.endswith(".0"):  # e.g. 95.0 -> 95
            number_string = number_string.split(".")[0]
        return res + [number_string + '\\% limits']

    def texValues(self, formatter, p, limit=2, refResults=None):
        if not isinstance(p, paramNames.paramInfo): param = self.parWithName(p)
        else: param = self.parWithName(p.name)
        if not param is None:
            lim = param.limits[limit - 1]
            sf = 3
            if lim.twotail:
                if not formatter.numberFormatter.plusMinusLimit(limit, lim.upper - param.mean, lim.lower - param.mean):
                    res, plus_str, _ = formatter.numberFormatter.namesigFigs(param.mean, param.err, param.err, wantSign=False)
                    res += r'\pm ' + plus_str
                else:
                    res, plus_str, minus_str = formatter.numberFormatter.namesigFigs(param.mean, lim.upper - param.mean, lim.lower - param.mean)
                    res += '^{' + plus_str + '}_{' + minus_str + '}'
            elif lim.onetail_upper:
                res = '< ' + formatter.numberFormatter.formatNumber(lim.upper, sf)
            elif lim.onetail_lower:
                res = '> ' + formatter.numberFormatter.formatNumber(lim.lower, sf)
            else: res = formatter.noConstraint
            if refResults is not None and res != formatter.noConstraint:
                refVal = refResults.parWithName(param.name)
                if refVal is not None:
                    delta = (param.mean - refVal.mean) / refVal.err
                    res += '\quad(%+.1f \\sigma)' % (delta)
            if self.hasBestFit:  # add best fit too
                rangew = (lim.upper - lim.lower) / 10
                bestfit = formatter.numberFormatter.namesigFigs(param.best_fit, rangew, -rangew)[0]
                return [res, bestfit]
            return [res]
        else: return None


class likeStats(paramResults):
    def loadFromFile(self, filename):
        textFileLines = self.fileList(filename)
        results = dict()
        for line in textFileLines:
            if len(line.strip()) == 0: break
            name, value = [x.strip() for x in line.split('=')]
            results[name] = float(value)
        self.logLike_sample = results.get('Best fit sample -log(Like)', None)
        self.logMeanInvLike = results.get('Ln(mean 1/like)', None)
        self.meanLogLike = results.get('mean(-Ln(like))', None)
        self.logMeanLike = results.get('-Ln(mean like)', None)


class convergeStats(paramResults):
    def loadFromFile(self, filename):
        textFileLines = self.fileList(filename)
        self.R_eigs = []
        for i in range(len(textFileLines)):
            if textFileLines[i].find('var(mean)') >= 0:
                for line in textFileLines[i + 1:]:
                    if len(line.strip()) == 0:break
                    try: self.R_eigs.append(line.split()[1])
                    except: self.R_eigs.append('1e30')
            elif 'Parameter auto-correlations' in textFileLines[i]:
                self.auto_correlation_steps = [int(s) for s in textFileLines[i + 2].split()]
                self.auto_correlations = []
                self.auto_correlation_pars = []
                for line in textFileLines[i + 3:]:
                    if len(line.strip()) == 0:break
                    items = line.split(None, len(self.auto_correlation_steps) + 1)
                    self.auto_correlation_pars.append(items[0])
                    self.auto_correlations.append([float(s) for s in items[1:-1]])


    def worstR(self):
        return self.R_eigs[len(self.R_eigs) - 1]


