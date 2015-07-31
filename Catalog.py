from __future__ import division

from ephem import Date
from collections import Mapping
import weakref
import numpy as np
import csv

class Catalog(object):
    def __init__(self, rows = None, orderedby = None, mode = 'ignore', verbose = False, **kwargs):
        '''
        Create a Catalog object from a generator or by reading in data from a csv file.
        
        If no value for rows is specified, this will create an empty catalog object.
        If creating a Catalog from points, rows should be an iterator over Point objects.
        If creating a Catalog from file, rows should be the path to the appropriate csv file.
        In this case, the columns to use are specified as keyword arguments, with the keyword
        being the name of the column and the value being a type to which to cast values in
        that column. Additional constants (values constant over all points) can be added as
        keyword arguments. The (optional) parameter orderedby specifies a column to use for 
        indexing and slicing. By default, this constructor ignores errors in typecasting,
        inserting the value None when an error is encountered. Setting verbose to True causes
        it to print a short message for each error encountered, and setting verbose to None will
        prevent it from catching errors. The mode parameter controls the treatment of unexpected
        attributes in new points. The default value 'ignore' will do nothing, neither adding 
        them to the catalog nor raising any exceptions. 'numeric' will try to cast such attributes
        to float, and 'raw' will add them in their current form.
        '''
        self.cast = {key: val for key, val in kwargs.iteritems() if callable(val)}
        self.const = {key: val for key, val in kwargs.iteritems() if not callable(val)}
        self.orderedby = orderedby
        self.functions = {}
        self.properties = {}
        self._points = []
        self._mode = mode
        if isinstance(rows, str):
            with open(rows) as csvcat:
                catdict = csv.DictReader(csvcat, dialect = 'excel')
                for row in catdict:
#                    print row
                    if 'date' in row.keys():
                        if row['date'][0] != '#': self._points.append(Point(self, mode = mode, verbose = verbose, **row))  # ignore comment lines
                    else:
                        self._points.append(Point(self, mode = mode, verbose = verbose, **row))
        elif rows is not None:
            for row in rows:
                if isinstance(row, Point):
                    if orderedby is None: self.orderedby = row.parents[-1]().orderedby
                    if len(self.cast) == 0: self.cast = row.parents[-1]().cast
                    row.parents.append(weakref.ref(self))
                    self._points.append(row)
                else:
                    if len(self.cast) == 0: self.cast = {key: type(row[key]) for key in row}
                    self._points.append(Point(self, mode = mode, verbose = verbose, **row))
        if orderedby is not None: self._points.sort(key=lambda pt: getattr(pt, orderedby))
    
    def __getitem__(self, key):
        if self.orderedby is None:
            if isinstance(key, int): return self._points[key]
            else: return Catalog(self._points[key])
        if isinstance(key, slice):
            start = self.cast[self.orderedby](key.start) if key.start is not None else None
            stop = self.cast[self.orderedby](key.stop) if key.stop is not None else None
            step = type(stop - start)(key.step) if key.step is not None else None
            if key.step is None: return Catalog(self._span(start, stop))
            else: return Catalog(pt for pt in self._span(start, stop) if (getattr(pt, self.orderedby) - start)/step % 1 == 0)
        else:
            for point in self._points:
                if getattr(point, self.orderedby) == self.cast[self.orderedby](key):
                    return point
               
    def __iter__(self):
        for point in self._points: yield point
    
    def __len__(self):
        return len(self._points)
    
    def __getattr__(self, name):
        if name in self.cast or name in self.properties or name in self.functions:
            return (getattr(point, name) for point in self._points)
        elif name.startswith('avg_'):
            truename = name[4:]
            if truename in self.cast or truename in self.properties:
                return np.mean([getattr(pt, truename) for pt in self._points])
            elif truename in self.functions:
                return lambda *args: np.mean([getattr(pt, truename)(*args) for pt in self._points])
            
    def __setattr__(self, name, val):
        if callable(val):
            if len(inspect.getargspec(val).args) == 1: self.add_property(self, name, val)
            else: self.add_function(self, name, val)
        else: object.__setattr__(self, name, val)
                    
    def __delattr__(self, name):
        if name in self.cast:
            for pt in self._points: delattr(pt, name)
            del self.cast[name]
        elif name in self.properties:
            for pt in self._points: delattr(pt, name)
            del self.properties[name]
        elif name in self.functions:
            for pt in self._points: delattr(pt, name)
            del self.functions[name]
        elif name in self.const:
            del self.const[name]
        else: object.__delattr__(self, name)
    def __del__(self):
        for point in self._points:
            for parent in point.parents: 
                if parent() is None: point.parents.remove(parent)
        
    def __add__(self, dct): self.add_point(**dct)
        
    def __getstate__(self): return self.__dict__
    
    def __setstate__(self, d): self.__dict__.update(d)
        
    def _addpt(self, point):
        point.parents.append(weakref.ref(self))
        self._points.append(point)
            
    def add_property(self, name, prop):
        '''
        Add a property to this Catalog.
        
        Properties are evaluated lazily and their values are stored for each groupby.
        The property should be specified as a function acting on a Point object. 
        After the property has been added, it will be accessible as an attribute of
        each point having the name specified.
        '''
        self.properties[name] = prop
        
    def add_function(self, name, func):
        '''
        Add a function to this Catalog.
        
        Functions are evaluated lazily and their values are stored for each point
        and combination of arguments. The first argument of the function specified
        should be a Point object, and all arguments should be hashable to allow lookup
        of argument tuples. After the property has been added, it will be accessible
        as a callable attribute of each point having the name specified and taking all
        arguments but the first.
        '''
        self.functions[name] = func
        
    def add_constant(self, name, const):
        '''
        Add a constant to this catalog.
        
        Constants are stored only once for the entire catalog.
        '''
        self.const[name] = const
        
    def add_point(self, verbose = False, **attrs):
        '''
        Add a point to this Catalog.
        
        Supply the attribute values for the point as keyword arguments.
        Set verbose to True to be warned if errors in casting occur, and to None
        to turn off this error-catching.
        '''
        self._points.append(Point(self, verbose = verbose, mode = self.mode, **attrs))
        
    def append(self, verbose = False, **attrs):
        '''
        Add a point to this Catalog.
        
        Supply the attribute values for the point as keyword arguments.
        Set verbose to True to be warned if errors in casting occur, and to None
        to turn off this error-catching.
        '''
        self.add_point(verbose, **attrs)
    
    def _span(self, start, stop):
        for point in self._points:
            current_index = getattr(point, self.orderedby)
            if start is not None and current_index < start: continue
            elif start is not None and current_index >= stop: break
            else: yield point
            
    def span(self, start, stop):
        '''
        Get a generator object that iterates over a range of points.
        
        Start and stop should be values of the attribute used for slicing, or integer
        indices if no such attribute was set in the initialization.
        '''
        if self.orderedby == None:
            return Catalog(self._points[start:stop])
        start = self.cast[self.orderedby](start)
        stop = self.cast[self.orderedby](stop)
        return self._span(start, stop)

    def successor(self, index):
        if self.orderedby is None: return self.__getitem__(index + 1)
        marker = False
        for point in self._points:
            if getattr(point, self.orderedby) == self.cast[self.orderedby](index): 
                marker = True
                continue
            if marker: return point
    
    def head(self, n):
        '''
        Get a generator object that iterates over the first n points in this Catalog.
        '''
        i = 0
        for point in self._points:
            if i >= n: break
            yield point
            i += 1
    
    def items(self, *args, **kwargs):
        '''
        Return a dictionary containing generators of attribute values.
        
        Only attributes specified in the arguments are included, except when
        there are no (non-keyword) arguments, in which case all attributes
        are included. Optionally, masks for the attribute names can be
        specified as keyword  arguments: for example, calling this function
        with dates = 'date' would name the list of date values 'dates' 
        rather than 'date'.
        '''
        if len(args) == 0:
            args = (name for name in set(self.cast) | set(self.properties) if not name.startswith('_'))
        mask = ((name, kwargs[name]) if name in kwargs else (name, name) for name in set(args) | set(kwargs))
        return {masked: (getattr(point, name) for point in self._points) for masked, name in mask}
    
    def itemlists(self, *args, **kwargs):
        '''
        Return a dictionary containing lists of attribute values.
        
        Only attributes specified in the arguments are included, except when
        there are no (non-keyword) arguments, in which case all attributes
        are included. Optionally, masks for the attribute names can be
        specified as keyword  arguments: for example, calling this function
        with dates = 'date' would name the list of date values 'dates' 
        rather than 'date'.
        '''
        if len(args) == 0:
            args = (name for name in set(self.cast) | set(self.properties) if not name.startswith('_'))
        mask = ((name, kwargs[name]) if name in kwargs else (name, name) for name in set(args) | set(kwargs))
        return {masked: [getattr(point, name) for point in self._points] for masked, name in mask}
    
    def mask(self, **kwargs):
        '''
        Return a dictionary containing generators of attribute values.
        
        Only attributes specified as keyword arguments are included.
        The keys in the returned dictionary are the specified keywords.
        '''
        return {masked: (getattr(point, name) for point in self._points) for masked, name in kwargs}
    
    def masklists(self, **kwargs):
        '''
        Return a dictionary containing lists of attribute values.
        
        Only attributes specified as keyword arguments are included.
        The keys in the returned dictionary are the specified keywords.
        '''
        return {masked: [getattr(point, name) for point in self._points] for masked, name in kwargs}

    def rename(self, *args, **kwargs):
        '''
        Change the name of an attribute of this catalog.
        
        To rename a single attribute, pass the name and then the
        new name as strings. To rename multiple attributes at once,
        pass a dictionary of newname: name pairs, or the equivalent
        as keyword arguments.
        '''
        if len(args) and isinstance(args[0], Mapping): names = args[0]
        elif len(kwargs): names = kwargs
        else: names = {args[1]: args[0]}
        for newname, name in names.iteritems():
            if name in self.cast:
                for pt in self._points:
                    setattr(pt, newname, getattr(pt, name))
                    delattr(pt, name)
                self.cast[newname] = self.cast[name]
                del self.cast[name]
            elif name in self.properties:
                for pt in self._points:
                    delattr(pt, name)
                self.properties[newname] = self.properties[name]
                del self.properties[name]
            elif name in self.functions:
                for pt in self._points:
                    delattr(pt, name)
                self.functions[newname] = self.functions[name]
                del self.functions[name]
            elif name in self.const:
                self.const[newname] = self.const[name]
                del self.const[name]
    
    def refactor(self, name, func):
        '''
        Apply func to the attribute name of each point in this catalog.
        '''
        if name in self.cast:
            for pt in self._points:
                setattr(pt, name, func(getattr(pt, name)))
            cast = self.cast[name]
            self.cast[name] = lambda x: func(cast(x))
        elif name in self.properties:
            for pt in self._points:
                setattr(pt, name, func(getattr(pt, name)))
            cast = self.properties[name]
            self.properties[name] = lambda x: func(cast(x))
        elif name in self.functions:
            for pt in self._points:
                delattr('_' + name + '_cache')
            cast = self.functions[name]
            self.functions[name] = lambda x: func(cast(x))
        elif name in self.const:
            self.const[name] = func(self.const[name])
        
    def groupby(self, key):
        '''
        Group the points in this Catalog in a specified manner.
        
        Groups points by their values of the attribute named by key, collecting
        points with the same value of the attribute into a catalog. Returns a
        dictionary of Catalogs indexed by attribute value. If key is instead
        specified as a slice object, the points with values of the orderedby
        attribute between start and stop will be grouped into bins of width step
        and indexed by their starting values.
        '''
        if isinstance(key, slice):
            step = self.cast[self.orderedby](key.step) if self.orderedby is not None else key.step
            groups = {}
            
            k = start
            for i, point in enumerate(self[start:stop]):
                ind = getattr(point, self.orderedby) if self.orderedby is not None else i
                if ind < start: break
                if ind >= k + step:
                    groups[k] = Catalog(orderedby = self.orderedby, **self.cast)
                    k = k + step
                groups[k]._addpt(point)
            return groups
        else:
            groups = {}
            
            for point in self._points:
                if callable(key): current_key = key(point)
                else: current_key = getattr(point, key)
                if not current_key in groups:
                    groups[current_key] = Catalog(orderedby = self.orderedby, **self.cast)
                groups[current_key]._addpt(point)
            return groups
    
    def orderby(self, attr):
        '''
        Change how this catalog is ordered.
        
        This affects dictionary-style keyed access.
        '''
        self.orderedby = attr
        if attr is not None: self._points.sort(key = lambda pt: getattr(pt, attr))
    
    def write(self, filename, fieldnames = None, writeprops = False, writeconsts = False, *names):
        '''
        Write this Catalog in csv format to the file at path filename.
        
        Specifying names in fieldnames or as additional arguments will cause only 
        attributes with those names to be written to the file, and in that specific
        order. (Otherwise order is arbitrary.) In the default mode, properties and
        constants will not be written to the file. To change this behavior, set 
        writeprops and/or writeconsts to True.
        '''
        if len(names) == 0: names = list(self.cast)
        if fieldnames is None: fieldnames = list(names)
        if writeprops: names += list(self.properties)
        if writeconsts: names += list(self.const)
        rows = [{name: getattr(pt, name) for name in fieldnames} for pt in self._points]
        with open(filename,'wb') as out_file:
            writer = csv.DictWriter(out_file, fieldnames, dialect='excel')
            writer.writerow(dict((fn,fn) for fn in fieldnames))
            for row in rows:
                writer.writerow(row)
                
    def writes(self, fieldnames = None, writeprops = False, writeconsts = False, *names):
        string = ''
        for pt in self._points:
            string += str(pt.date)+','+str(pt.ra)+','+str(pt.dec)+','+str(pt.expnum)+','+str(pt.exptime)+','+pt.band+','+str(pt.ccdnum)+','+\
            str(pt.mag)+','+str(pt.objid)+'\n'
        return string
    
class Point(object):
    def __init__(self, parent, mode = None, verbose = False, **attributes):
        self.parents = [weakref.ref(parent)]
        for name, val in attributes.iteritems():
            sname = name.strip()
            try: 
                setattr(self, sname, None if val is None else (val if parent.cast[sname] == type(None) else parent.cast[sname](val)))
            except KeyError:
                if mode == 'ignore': pass
                elif mode == 'raw':
                    parent.cast.update({sname: type(val)})
                    setattr(self, sname, val)
                elif mode == 'numeric':
                    parent.cast.update({sname: float})

    @property
    def attributes(self):
        attributes = set()
        for parent in self.parents:
            attributes = attributes | set(name for name in parent().cast)
        for name in attributes:
            yield name
    
    @property
    def properties(self):
        for parent in self.parents:
            for name, prop in parent().properties.iteritems(): yield name, prop
    
    @property
    def functions(self):
        for parent in self.parents:
            for name, func in parent().functions.iteritems(): yield name, func
                
    @property
    def constants(self):
        for parent in self.parents:
            for name, const in parent().const.iteritems(): yield name, const

    def __iter__(self):
        for name in self.attributes: yield name
        for name, prop in self.properties: yield name
        for name, func in self.functions: yield name
        for name, const in self.constants: yield name
        
    def __getattr__(self, name):
        for key, prop in self.properties:
            if key == name:
                rval = prop(self)
                setattr(self, name, rval)
                return rval
        for key, func in self.functions:
            if key == name:
                setattr(self, '_' + name + '_cache', {})
                def wrapper(*args):
                    cache = getattr(self, '_' + name + '_cache')
                    if args in cache:
                        return cache[args]
                    else:
                        rval = func(self, *args)
                        cache[args] = rval
                        return rval
                setattr(self, name, wrapper)
                return wrapper
        for key, const in self.constants:
            if key == name: return const
    
    def __delattr__(self, name):
        for key, func in self.functions:
            if key == name:
                object.__delattr__(self, '_' + name + '_cache')
                break
        object.__delattr__(self, name)
    
    def __repr__(self):
        return "Catalog.Point(" + ', '.join(name + '=' + repr(getattr(self, name)) for name in self) + ')'

class TimeDelta(float):
    def __new__(cls, arg):
        if isinstance(arg, str):
            if len(arg.split()) == 1:
                return TimeDelta(float(arg))
            delta = 0
            for string in arg.split(','):
                val, unit = string.split()
                if unit in ['Y','a','y','yr','year','years']: unit = 365.25
                elif unit in ['W','w','week','weeks']: unit = 7
                elif unit in ['D','d','day','days']: unit = 1
                elif unit in ['h','hr','hour','hours']: unit = 1/24
                elif unit in ['m','min','minute','minutes']: unit = 1/1440
                elif unit in ['s','sec','second','seconds']: unit = 1/86400
                elif unit in ['ms','millisecond','milliseconds']: unit = 1/86400000
                elif unit in ['us','microsecond','microseconds']: unit = 1/86400000000
                delta += float(val) * unit
            return float.__new__(cls, delta)
        elif isinstance(arg, float):
            return float.__new__(cls, arg)
        elif isinstance(arg, np.timedelta64):
            return TimeDelta(str(arg))
    
    def __str__(self):
        formatter = '{0:.4g}'
        if self < TimeDelta('1 ms'): return formatter.format(86400000000*self) + ' us'
        elif self < TimeDelta('1 s'): return formatter.format(86400000*self) + ' ms'
        elif self < TimeDelta('1 m'): return formatter.format(86400*self) + ' sec'
        elif self < TimeDelta('1 h'): return str(int(1440*self % 60)) + ' min, ' + formatter.format(86400*self % 60) + ' sec'
        elif self < TimeDelta('1 D'): return str(int(24*self % 24)) + ' hr, ' + formatter.format(1440*self % 60) + ' min'
        elif self < TimeDelta('1 Y'): return str(int(self)) + ' d, ' + formatter.format(24*self % 24) + ' hr'
        elif self < TimeDelta('10 Y'): return str(int(self / 365.25)) + ' yr, ' + str(int(self % 365.25)) + ' d'
        else: return str(int(self / 365.25)) + ' yr'
    
    def __sub__(self, other): return TimeDelta(super(TimeDelta, self).__sub__(other))
    
    def __add__(self, other): return TimeDelta(super(TimeDelta, self).__add__(other))
    
    def __mult__(self, other): return TimeDelta(super(TimeDelta, self).__mult__(other))

    @property
    def day(self): return int(self)

    @property
    def hour(self): return int(24 * (self % 1))

    @property
    def minute(self): return int(1440 * (self % (1/24)))

    @property
    def second(self): return int(86400 * (self % (1/1440)))

    @property
    def millisecond(self): return 86400000 * (self % (1/86400))

class DateTime(Date):
    def __new__(cls, *args): return Date.__new__(cls, *args)
    
    def __sub__(self, other): return TimeDelta(super(DateTime, self).__sub__(other))
    
    def __add__(self, other): return DateTime(super(DateTime, self).__add__(other))

    @property
    def day(self): return int(self)

    @property
    def hour(self): return int(24 * (self % 1))

    @property
    def minute(self): return int(1440 * (self % (1/24)))

    @property
    def second(self): return int(86400 * (self % (1/1440)))

    @property
    def millisecond(self): return 86400000 * (self % (1/86400))
            

