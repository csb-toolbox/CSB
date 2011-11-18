"""
Plotting facilities, based on Python's MPL library.

The L{Chart} object is a facade which provides intuitive access to MPL's plotting
objects. The following examples shows a typical use of L{Chart}:

    1. Instantiation
        
    >>> chart = Chart()
    
    2. Adding plots (equivalent to MPL's subplots)
    
    >>> chart.plots.add(1, 2, 1)
    Plot (matplotlib.Axes)
    >>> chart.plots.add(1, 2, 2)
    Plot (matplotlib.Axes)
    
    3. Plotting
    
    >>> chart.plots[0].hist(...)
    >>> chart.plots[0].set_title('T')
        
    >>> chart.plots[1].scatter(...)
    >>> chart.plots[1].set_xlabel('X')
    
    4. Using the GUI
    
    >>> chart.show()
    >>> chart.hide()
    
    5. Saving as image
    
    >>> chart.save(filename)
    
If using the GUI, do not forget to dispose the chart at the end:

    >>> chart.dispose()
    
or simply use the chart in a context manager:

    >>> with Chart() as chart:
            chart...
"""

import time

from abc import ABCMeta, abstractmethod

from Queue import Queue
from threading import Thread, Event

from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg


class Backend(Thread):
    """
    Abstract class defining the behavior of all Chart GUI backends.
    
    Each backend is a 'daemon' that runs in a new thread and handles
    all GUI requests from any L{Chart} instance. A backend service must
    behave as a singleton - only one service of a given kind may exist at a
    given point in time. L{Chart} clients will request GUI operations
    on specific figures, the Backend therefore must keep track of all
    windows opened, as well as the figure-to-window mapping.
    """
    
    __metaclass__ = ABCMeta
    _instances = {}
    
    @staticmethod
    def get(backend, started=True):
        """
        Backend factory, ensures one instance per subclass. 
        
        @param backend: one of the L{Backend} subclasses
        @type backend: type
        @param started: if True, ensure that the service is running
        @type started: bool
        
        @return: an instance of the backend. The returned service
                 instance may need to be started.
        @rtype: L{Backend}
        """

        if not issubclass(backend, Backend):
            raise TypeError(backend)
        
        if backend in Backend._instances:
            instance = Backend._instances[backend]
        else:
            instance = backend()

        if started and not instance.started:
            instance.start()  
        return instance
    
    def __init__(self):
        
        name = self.__class__
        if name in Backend._instances:
            raise RuntimeError('Backend {0} has already been initialized'.format(name))
        else:
            Backend._instances[name] = self         
        
        super(Backend, self).__init__()
        
        self._figures = {}
        self._queue = Queue()
        self._started = Event()        
        self._running = Event()
                
        self.setDaemon(True)
        
    @property
    def started(self):
        return self._started.isSet()
    
    @property
    def running(self):
        return self._running.isSet()    
                
    @abstractmethod
    def _initapp(self):
        """
        Create an instance of the GUI application.
        """
        pass
    
    @abstractmethod
    def _mainloop(self):
        """
        Enter the GUI main loop.
        """
        pass
    
    @abstractmethod
    def _exit(self):
        """
        Delete all frames, exit the GUI main loop and perform any cleanup
        needed in order to unblock the thread that started the main loop.
        """
        pass    
    
    @abstractmethod
    def _add(self, figure):
        """
        Handle an 'Add new figure' event
        """
        pass
    
    @abstractmethod    
    def _show(self, figure):
        """
        Handle a 'Show existing figure' event
        """
        pass

    @abstractmethod    
    def _hide(self, figure):
        """
        Handle a 'Hide existing figure' event
        """
        pass

    @abstractmethod        
    def _destroy(self, figure):
        """
        Handle a 'Delete existing figure' event
        """        
        pass

    @abstractmethod
    def _invoke(self, callable, *args):
        """
        Pass a GUI message: invoke C{callable} in a thread-safe way
        """
        pass
    
    def invoke(self, callable, *args):
        """
        Invoke an asynchronous GUI operation (in a thread-safe way) 
        """
        if not self._running.isSet():
            raise RuntimeError('The backend service is not running')
        else:
            self._invoke(callable, *args)
                
    def add(self, figure):
        self.invoke(self._add, figure)

    def show(self, figure):
        self.invoke(self._show, figure)
    
    def hide(self, figure):
        self.invoke(self._hide, figure)
        
    def destroy(self, figure, wait=False):
        
        has_figure = (figure in self._figures)
        self.invoke(self._destroy, figure)
        
        if has_figure and wait:
            while figure in self._figures:
                pass       
                
    def start(self):
        """
        Start the Backend service. This method can be called only once.
        """
        try:
            super(Backend, self).start()
        
            while not self._running.isSet():
                time.sleep(0.05)
                
        except BaseException:
            raise RuntimeError("Failed to start the backend service")
    
    def run(self):
        """
        Main service method, automatically called by C{start}.
        """
        self._started.set()
        
        self._initapp()
        self._running.set()
        
        self._mainloop()

        self._running.clear()
        self._started.clear()
        
    def stop(self):
        """
        Stop the Backend service. The Backend object can be safely
        disposed afterwards.
        """
        self._exit()
        self._figures = {}
        self.join()
        self._running.clear()
        self._started.clear()
        del Backend._instances[self.__class__]
    
    def client_disposed(self, client):
        """
        Fired when a client is being deleted. Will stop the service if no
        active clients are remaining.
        """
        if self._figures is None or len(self._figures) == 0:
            self.stop()
            
    def __del__(self):
        if self._started.isSet():
            self.stop()
    
class WxBackendImpl(Backend):
    """
    WxPython L{Backend} implementor. 
    """
    
    _wxapp = None
    
    def __init__(self):
        
        import wx
        from matplotlib.backends.backend_wx import FigureFrameWx
        
        self.wx = wx
        self.FigureFrameWx = FigureFrameWx
        
        super(WxBackendImpl, self).__init__()
    
    @property
    def _app(self):
        if WxBackendImpl._wxapp is None:
            WxBackendImpl._wxapp = self.wx.PySimpleApp()
        return WxBackendImpl._wxapp
                    
    def _initapp(self):
      
        dummy = self._app
        frame = self.wx.Frame(None)
        frame.Show()
        frame.Hide()
        
    def _mainloop(self):
        
        self._app.MainLoop()
        
    def _add(self, figure):
                
        wx = self.wx
        FigureFrameWx = self.FigureFrameWx
        
        if figure not in self._figures:
            
            frame = FigureFrameWx(figure._figure_number, figure)
            frame.Show()
            frame.Bind(wx.EVT_ACTIVATE, lambda e: e.GetEventObject().Layout())
            frame.Bind(wx.EVT_CLOSE, lambda e: self.invoke(self._destroy, figure))

            self._figures[figure] = frame

    def _show(self, figure):
        
        if figure not in self._figures:
            self._add(figure)
        
        self._figures[figure].Show()
            
    def _hide(self, figure):
        
        if figure in self._figures:
            self._figures[figure].Hide()                

    def _destroy(self, figure):
        
        if figure in self._figures:
            frame = self._figures[figure]
            if not frame.IsBeingDeleted():
                frame.Destroy()            
            del self._figures[figure]
                
    def _invoke(self, callable, *args):
        
        wx = self.wx
        wx.CallAfter(callable, *args)
                    
    def _exit(self):
        
        for frame in self._figures.values():
            if not frame.IsBeingDeleted():
                frame.Destroy()
        self._app.Exit()
        
class Backends(object):
    """
    Enumeration of chart backends.
    """
    
    WX_WIDGETS = WxBackendImpl
        

class PlotsCollection(object):
    """
    A list-like collection of all plots in the chart (0-based).
    """
    
    def __init__(self, figure):
        
        self._plots = []
        self._figure = figure
        
    def add(self, rows=1, columns=1, cell=1):
        """
        Add a new plot to the figure. By default the new plot will span the
        whole figure (window). The optional parameters can be used to position
        the plot explicitly by dividing the figure virtually into a grid and
        specifying the location of the plot in that grid. All cells in the
        grid are numbered sequentially, left-to-right, starting from 1.
        
        @param rows: number of "rows" in the grid
        @type rows: int
        @param columns: number of "columns" in the grid
        @type columns: int
        @param cell: in which "cell" to put the new plot
        @type cell: int
        """
        for i in [rows, columns, cell]:
            if i < 1:
                raise ValueError(i)
            
        if not cell <= rows * columns:
            raise ValueError('Position {0} is greater than the number'
                             ' of cells in a {0}x{1} grid'.format(rows, columns, cell))
            
        plot = self._figure.add_subplot(rows, columns, cell)
        self._plots.append(plot)
        
        return plot
    
    def __getitem__(self, i):
        try:
            return self._plots[i]
        except IndexError:
            raise IndexError("No such plot number: {0}".format(i))
        
    def __len__(self):
        return len(self._plots)
    
    def __iter__(self):
        return iter(self._plots)
    
              
class Chart(object):
    """
    Simple and clean facade to Matplotlib's plotting API.
    
    A chart instance abstracts a plotting device, on which one or
    multiple related plots can be drawn. Charts can be exported as images, or
    visualized interactively. Each chart instance will always open in its own
    GUI window, and this window will never block the execution of the rest of
    the program, or interfere with other L{Chart}s.
    The GUI can be safely opened in the background and closed infinite number
    of times, as long as the client program is still running.
    
    @param number: chart number; by default this a L{Chart.AUTONUMBER}
    @type number: int or None
    @param title: chart master title
    @type title: str
    
    @note: additional arguments are passed directly to Matplotlib's Figure
           constructor. 
    """

    AUTONUMBER = None
    
    _serial = 0
    
    
    def __init__(self, number=None, title='', backend=Backends.WX_WIDGETS, *fa, **fk):
        
        if number == Chart.AUTONUMBER:
            Chart._serial += 1
            number = Chart._serial
            
        self._number = int(number)
        self._title = str(title)
        self._figure = Figure(*fa, **fk)
        self._figure._figure_number = self._number
        self._figure.suptitle(self._title)
        self._beclass = backend
        self._hasgui = False
        self._canvas = FigureCanvasAgg(self._figure)
        self._plots = PlotsCollection(self._figure)
    
    def __getitem__(self, i):
        if i in self._plots:
            return self._plots[i]
        else:
            raise KeyError('No such plot number: {0}'.format(i))
        
    def __enter__(self):
        return self
    
    def __exit__(self, *a, **k):
        self.dispose()
    
    @property
    def _backend(self):
        return Backend.get(self._beclass, started=True)
  
    @property
    def title(self):
        return self._title
        
    @property
    def number(self):
        return self._number
    
    @property
    def plots(self):
        return self._plots
    
    @property     
    def figure(self):
        return self._figure
            
    def show(self):
        """
        Show the GUI window (non-blocking).
        """
        if not self._hasgui:
            self._backend.add(self._figure)
            self._hasgui = True
            
        self._backend.show(self._figure)
                
    def hide(self):
        """
        Hide (but do not dispose) the GUI window.
        """
        self._backend.hide(self._figure)
        
    def dispose(self):
        """
        Dispose the GUI interface. Must be called at the end if any
        chart.show() calls have been made. Automatically called if using
        the chart in context manager ("with" statement).
        
        @note: Failing to call this method if show() has been called at least
        once may cause backend-related errors.
        """
        service = self._backend
        
        if service and service.running:
            service.destroy(self._figure, wait=True)
            service.client_disposed(self)    
        
    def save(self, file, format='png', dpi=None, *a, **k):
        """
        Save all plots to an image.
        """
        self._canvas.print_figure(file, format=format, dpi=dpi, *a, **k)
