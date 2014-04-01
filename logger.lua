local class   = require 'class'
local colors  = require 'ansicolors'

local CommandLineLogger = class.class 'CommandLineLogger'

function CommandLineLogger:__init__(class_name)
   self._class_name = class_name
   self._enable_colors = true
end
function CommandLineLogger:log_message(funcname, msg, indent)
   local sep = self._class_name and ':' or '->'
   if self._enable_colors then
      print(colors(string.rep('  ', indent or 0)..
		   '%{magenta}[%{green}'..(self._class_name or '')..
		   '%{cyan}'..sep..'%{underline}'..funcname..
		   '%{reset}%{magenta}]%{black} '..msg))
   else
      print(string.rep('  ', indent or 0)..
	    '['..(self._class_name or '')..sep..funcname..'] '..msg)
   end
end
function CommandLineLogger:disable_colors()
   self._enable_colors = false
end
function CommandLineLogger:enable_colors()
   self._enable_colors = true
end

return {CommandLineLogger=CommandLineLogger}

