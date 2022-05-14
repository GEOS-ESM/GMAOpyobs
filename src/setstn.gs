function setstn(arg)

   if (arg='')
     say ''
     say 'setstn - set ensemble environment and return station name'
     say 'Usage:'
     say '        setstn e'
     say 'where e is the ensemble number'
     say ''
     return
   else
     e = arg
     'set e ' e
      if (e < 10)
         aname = 'Station_00'e
      else; if (e<100)
         aname = 'Station_0'e
      else; if (e<1000)
         aname = 'Station_'e
      else
         say 'too many stations'
         return
      endif; endif; endif
   endif

*  Find station name
*  -----------------
   'q attr'
   attrs = result
   i = 3
   line = sublin(attrs,i)
   while ( line != '' )
     if ( subwrd(line,1)='global' )
        if ( subwrd(line,2)='String' )
          name = subwrd(line,3)
          if ( name=aname )
               stn = subwrd(line,4)' 'subwrd(line,5)' 'subwrd(line,6)
               return stn
          endif
        endif
     endif
     i = i + 1
     line = sublin(attrs,i)

   endwhile

return '***error***'

