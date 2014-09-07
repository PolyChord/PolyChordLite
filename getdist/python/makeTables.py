import os, batchJobArgs, ResultObjs, paramNames, planckStyle, copy


Opts = batchJobArgs.batchArgs('Make pdf tables from latex generated from getdist outputs', importance=True, converge=True)
Opts.parser.add_argument('latex_filename')
Opts.parser.add_argument('--limit', type=int, default=2)
Opts.parser.add_argument('--all_limits', action='store_true')

Opts.parser.add_argument('--bestfitonly', action='store_true')
Opts.parser.add_argument('--nobestfit', action='store_true')
Opts.parser.add_argument('--no_delta_chisq', action='store_true')
Opts.parser.add_argument('--delta_chisq_paramtag', default=None)
Opts.parser.add_argument('--changes_from_datatag', default=None)
Opts.parser.add_argument('--changes_from_paramtag', default=None)
Opts.parser.add_argument('--changes_adding_data', default=None)

# this is just for the latex labelsm set None to use those in chain .paramnames
Opts.parser.add_argument('--paramNameFile', default='clik_latex.paramnames')
Opts.parser.add_argument('--paramList', default=None)
Opts.parser.add_argument('--blockEndParams', default=None)
Opts.parser.add_argument('--columns', type=int, nargs=1, default=3)
Opts.parser.add_argument('--compare', nargs='+', default=None)

Opts.parser.add_argument('--titles', default=None)  # for compare plots
Opts.parser.add_argument('--forpaper', action='store_true')
Opts.parser.add_argument('--separate_tex', action='store_true')
Opts.parser.add_argument('--header_tex', default=None)
Opts.parser.add_argument('--height', default="9in")
Opts.parser.add_argument('--width', default="11in")

(batch, args) = Opts.parseForBatch()

if args.blockEndParams is not None: args.blockEndParams = args.blockEndParams.split(';')

if args.paramList is not None: args.paramList = paramNames.paramNames(args.paramList)

if args.forpaper: formatter = planckStyle.planckStyleTableFormatter()
else: formatter = None

def texEscapeText(string):
    return string.replace('_', '{\\textunderscore}')

def getTableLines(content, referenceDataJobItem=None):
    if referenceDataJobItem is not None: refResults = referenceDataJobItem.result_marge
    else: refResults = None
    return ResultObjs.resultTable(args.columns, [content], blockEndParams=args.blockEndParams,
                         formatter=formatter, paramList=args.paramList, limit=args.limit, refResults=refResults).lines

def paramResultTable(jobItem, deltaChisqJobItem=None, referenceDataJobItem=None):
    if deltaChisqJobItem is not None and deltaChisqJobItem.name == jobItem.name: deltaChisqJobItem = None
    if referenceDataJobItem is not None:
        if (args.changes_from_paramtag is None and referenceDataJobItem.normed_data == jobItem.normed_data
         or args.changes_from_paramtag is not None and referenceDataJobItem.name == jobItem.name):
            referenceDataJobItem = None
    tableLines = []
    caption = []
    jobItem.loadJobItemResults(paramNameFile=args.paramNameFile, bestfit=not args.nobestfit, bestfitonly=args.bestfitonly)
    bf = jobItem.result_bestfit
    if not bf is None:
        caption.append(' Best-fit $\\chi^2_{\\rm eff} = ' + ('%.2f' % (bf.logLike * 2)) + '$')
        if deltaChisqJobItem is not None:
            bf_ref = deltaChisqJobItem.result_bestfit
            if bf_ref is not None: caption.append('$\\Delta \\chi^2_{\\rm eff} = ' + ('%.2f' % ((bf.logLike - bf_ref.logLike) * 2)) + '$')

    if args.bestfitonly:
        if bf is not None: tableLines += getTableLines(bf)
    else:
        likeMarge = jobItem.result_likemarge
        if likeMarge is not None and likeMarge.meanLogLike is not None:
                caption.append('$\\bar{\\chi}^2_{\\rm eff} = ' + ('%.2f' % (likeMarge.meanLogLike * 2)) + '$')
                if deltaChisqJobItem is not None:
                    likeMarge_ref = deltaChisqJobItem.result_likemarge
                    if likeMarge_ref is not None and likeMarge_ref.meanLogLike is not None:
                        delta = likeMarge.meanLogLike - likeMarge_ref.meanLogLike
                        caption.append('$\\Delta\\bar{\\chi}^2_{\\rm eff} = ' + ('%.2f' % (delta * 2)) + '$')
        if jobItem.result_converge is not None: caption.append('$R-1 =' + jobItem.result_converge.worstR() + '$')
        if jobItem.result_marge is not None: tableLines += getTableLines(jobItem.result_marge, referenceDataJobItem)
    tableLines.append('')
    if not args.forpaper: tableLines.append("; ".join(caption))
    if not bf is None and not args.forpaper:
        tableLines.append('')
        tableLines.append('$\\chi^2_{\\rm eff}$:')
        if deltaChisqJobItem is not None: compChiSq = deltaChisqJobItem.result_bestfit
        else: compChiSq = None
        for kind, vals in bf.sortedChiSquareds():
            tableLines.append(kind + ' - ')
            for (name, chisq) in vals:
                line = '  ' + texEscapeText(name) + ': ' + ('%.2f' % chisq) + ' '
                if compChiSq is not None:
                    comp = compChiSq.chiSquareForKindName(kind, name)
                    if comp is not None: line += '($\Delta$ ' + ('%.2f' % (chisq - comp)) + ') '
                tableLines.append(line)
    return tableLines

def compareTable(jobItems, titles=None):
    for jobItem in jobItems:
        jobItem.loadJobItemResults(paramNameFile=args.paramNameFile, bestfit=not args.nobestfit, bestfitonly=args.bestfitonly)
        print jobItem.name
    if titles is None: titles = [jobItem.datatag for jobItem in jobItems if jobItem.result_marge is not None]
    else: titles = titles.split(';')
    return ResultObjs.resultTable(1, [jobItem.result_marge for jobItem in jobItems if jobItem.result_marge is not None],
               formatter=formatter, limit=args.limit, titles=titles, blockEndParams=args.blockEndParams, paramList=args.paramList).lines


items = Opts.sortedParamtagDict(chainExist=not args.bestfitonly)

if args.all_limits:
    limits = [1, 2, 3]
else: limits = [args.limit]


if args.changes_from_paramtag is not None:
    if args.changes_from_datatag is not None:
        raise Exception('You cannot have both changes_from_paramtag and changes_from_datatag')
    if args.delta_chisq_paramtag is not None and args.delta_chisq_paramtag != args.changes_from_paramtag:
        raise Exception('when using changes_from_paramtag, delta_chisq_paramtag is set equal to that')
    if args.no_delta_chisq:
        raise Exception('when using changes_from_paramtag cannot have no_delta_chisq')
    args.delta_chisq_paramtag = args.changes_from_paramtag

baseJobItems = dict()
for paramtag, parambatch in items:
    isBase = len(parambatch[0].param_set) == 0
    for jobItem in parambatch:
        if (args.delta_chisq_paramtag is None and
            isBase and not args.no_delta_chisq  or args.delta_chisq_paramtag is not None and jobItem.paramtag == args.delta_chisq_paramtag):
                referenceJobItem = copy.deepcopy(jobItem)
                referenceJobItem.loadJobItemResults(paramNameFile=args.paramNameFile)
                baseJobItems[jobItem.normed_data] = referenceJobItem


for limit in limits:
    args.limit = limit

    outfile = args.latex_filename
    if args.all_limits: outfile += '_limit' + str(limit)
    if outfile[-4:] != '.tex': outfile += '.tex'

    lines = []
    if not args.forpaper:
        lines.append('\\documentclass[10pt]{article}')
        lines.append('\\usepackage{fullpage}')
        lines.append('\\usepackage[pdftex]{hyperref}')
        lines.append('\\usepackage[paperheight=' + args.height + ',paperwidth=' + args.width + ',margin=0.8in]{geometry}')
        lines.append('\\renewcommand{\\arraystretch}{1.5}')
        lines.append('\\begin{document}')
        if args.header_tex is not None:
            lines.append(open(args.header_tex, 'r').read())
        lines.append('\\tableofcontents')

    # set of baseline results, e.g. for Delta chi^2

    for paramtag, parambatch in items:
        isBase = len(parambatch[0].param_set) == 0
        if not args.forpaper:
            if isBase: paramText = 'Baseline model'
            else: paramText = texEscapeText("+".join(parambatch[0].param_set))
            section = '\\newpage\\section{ ' + paramText + '}'
        else: section = ''
        if args.compare is not None:
            compares = Opts.filterForDataCompare(parambatch, args.compare)
            if len(compares) == len(args.compare):
                lines.append(section)
                lines += compareTable(compares, args.titles)
            else: print 'no matches for compare: ' + paramtag
        else:
            lines.append(section)
            theseItems = [jobItem for jobItem in parambatch
                if (os.path.exists(jobItem.distPath) or args.bestfitonly) and (args.converge == 0 or jobItem.hasConvergeBetterThan(args.converge))]

            referenceDataJobItem = None
            if args.changes_from_datatag is not None:
                for jobItem in theseItems:
                    if jobItem.normed_data == args.changes_from_datatag or jobItem.datatag == args.changes_from_datatag:
                        referenceDataJobItem = copy.deepcopy(jobItem)
                        referenceDataJobItem.loadJobItemResults(paramNameFile=args.paramNameFile, bestfit=args.bestfitonly)
            if args.changes_adding_data is not None:
                baseJobItems = dict()
                refItems = []
                for jobItem in theseItems:
                    if jobItem.data_set.hasName(args.changes_adding_data):
                        jobItem.normed_without = "_".join(sorted([x for x in jobItem.data_set.names if not x == args.changes_adding_data]))
                        refItems.append(jobItem.normed_without)
                    else: jobItem.normed_without = None
                for jobItem in theseItems:
                    if jobItem.normed_data in refItems:
                        referenceJobItem = copy.deepcopy(jobItem)
                        referenceJobItem.loadJobItemResults(paramNameFile=args.paramNameFile, bestfit=args.bestfitonly)
                        baseJobItems[jobItem.normed_data] = referenceJobItem

            for jobItem in theseItems:
                    if not args.forpaper: lines.append('\\subsection{ ' + texEscapeText(jobItem.name) + '}')
                    if args.changes_adding_data is not None:
                        if jobItem.normed_without is not None:
                            referenceDataJobItem = baseJobItems.get(jobItem.normed_without, None)
                        else: referenceDataJobItem = None
                        referenceJobItem = referenceDataJobItem
                    else: referenceJobItem = baseJobItems.get(jobItem.normed_data, None)
                    if args.changes_from_paramtag is not None:
                        referenceDataJobItem = referenceJobItem

                    tableLines = paramResultTable(jobItem, referenceJobItem, referenceDataJobItem)
                    if args.separate_tex: ResultObjs.textFile(tableLines).write(jobItem.distRoot + '.tex')
                    lines += tableLines

    if not args.forpaper: lines.append('\\end{document}')

    (outdir, outname) = os.path.split(outfile)
    if len(outdir) > 0 and not os.path.exists(outdir): os.makedirs(outdir)
    ResultObjs.textFile(lines).write(outfile)
    (root, _) = os.path.splitext(outfile)

    if not args.forpaper:
        print 'Now converting to PDF...'
        delext = ['aux', 'log', 'out', 'toc']
        if len(outdir) > 0: dodir = 'cd ' + outdir + '; '
        else: dodir = '';
        os.system(dodir + 'pdflatex ' + outname)
        # #again to get table of contents
        os.system(dodir + 'pdflatex ' + outname)
        # and again to get page numbers
        os.system(dodir + 'pdflatex ' + outname)
        for ext in delext:
            if os.path.exists(root + '.' + ext):
                os.remove(root + '.' + ext)


