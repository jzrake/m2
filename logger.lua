local class   = require 'class'
local colors  = require 'ansicolors'

local CommandLineLogger = class.class 'CommandLineLogger'
CommandLineLogger._rc = {
   enable_colors=true,
   message_body_color='white'
}
function CommandLineLogger:__init__(class_name)
   self._class_name = class_name
   self._enable_colors = CommandLineLogger._rc.enable_colors
   self._message_body_color = CommandLineLogger._rc.message_body_color
end
function CommandLineLogger:log_message(funcname, msg, indent)
   local sep = self._class_name and ':' or '->'
   if self._enable_colors then
      print(colors(string.rep('  ', indent or 0)..
		   '%{magenta}[%{green}'..(self._class_name or '')..
		   '%{cyan}'..sep..'%{underline}'..funcname..
		   '%{reset}%{magenta}]%{'..
		   self._message_body_color..'} '..msg))
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
