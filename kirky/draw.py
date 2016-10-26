from pyx import path, deco, text, color
from copy import copy


def DrawEdge(edge, canvas):
    head = [float(element) for element in edge.head]
    tail = [float(element) for element in edge.tail]
    reversed = False
    if tail[0] - head[0] > 0:
        reversed = True
    for i in range(0, len(head)):
        head[i] *= 4
        tail[i] *= 4
    weight = float(edge.weight)
    if weight == 0:
        return
    string = '%s' % weight
    """
    numerator = weight.numerator
    #if numerator == 0:
    #    return
    denominator = weight.denominator
    string = '%s/%s' % (numerator, denominator)
    """
    if not reversed:
        if string and string[0] == '-':
            canvas.stroke(path.line(tail[0], tail[1], head[0], head[1]), [deco.earrow, deco.curvedtext(string, textattrs=[text.halign.left, text.vshift.mathaxis, text.size.tiny], arclenfrombegin=1.2, exclude=0.1), color.rgb.red])
        elif string and string[0] != '-':
            canvas.stroke(path.line(tail[0], tail[1], head[0], head[1]), [deco.earrow, deco.curvedtext(string, textattrs=[text.halign.left, text.vshift.mathaxis, text.size.tiny], arclenfrombegin=1.2, exclude=0.1)])
    else:
        if string and string[0] == '-':
            canvas.stroke(path.line(head[0], head[1], tail[0], tail[1]), [deco.barrow, deco.curvedtext(string, textattrs=[text.halign.left, text.vshift.mathaxis, text.size.tiny], arclenfrombegin=1.2, exclude=0.1), color.rgb.red])
        elif string and string[0] != '-':
            canvas.stroke(path.line(head[0], head[1], tail[0], tail[1]), [deco.barrow, deco.curvedtext(string, textattrs=[text.halign.left, text.vshift.mathaxis, text.size.tiny], arclenfrombegin=1.2, exclude=0.1)])
