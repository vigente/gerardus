/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file  spxlpfread.cpp
 * @brief Read LP format files.
 */
// #define DEBUGGING 1

#include <assert.h>
#include <stdio.h>
#include <ctype.h>
#include <iostream>

#include "spxdefines.h"
#include "spxlp.h"
#include "spxout.h"
#include "exceptions.h"

/* The manual says the maximum allowed line length is 255 characters,
 * but CPLEX does not complain if the lines are longer.
 */
#define MAX_LINE_LEN  8192     ///< maximum length of a line (8190 + \\n + \\0)

namespace soplex
{
/// Is \p c a \c space, \c tab, \c nl or \c cr ?
static inline bool isSpace(int c)
{
   return (c == ' ') || (c == '\t') || (c == '\n') || (c == '\r'); 
}

/// Is there a number at the beginning of \p s ?
static bool isValue(const char* s)
{
   return ((*s >= '0') && (*s <= '9'))  
      || (*s == '+') || (*s == '-') || (*s == '.');
}

/// Is there a possible column name at the beginning of \p s ?
static bool isColName(const char* s)
{
   // strchr() gives a true for the null char.
   if (*s == '\0')
      return false;

   return ((*s >= 'A') && (*s <= 'Z'))
      || ((*s >= 'a') && (*s <= 'z'))
      || (strchr("!\"#$%&()/,;?@_'`{}|~", *s) != 0);
}

/// Is there a comparison operator at the beginning of \p s ?
static bool isSense(const char* s)
{
   return (*s == '<') || (*s == '>') || (*s == '=');
}

static bool isInfinity(const char* s)
{
   return ((s[0] == '-') || (s[0] == '+'))
      && (tolower(s[1]) == 'i') 
      && (tolower(s[2]) == 'n') 
      && (tolower(s[3]) == 'f');
}

static bool isFree(const char* s)
{
   return (tolower(s[0]) == 'f') 
      && ( tolower(s[1]) == 'r') 
      && ( tolower(s[2]) == 'e') 
      && ( tolower(s[3]) == 'e');
}

/// Read the next number and advance \p pos.
/** If only a sign is encountered, the number is assumed to 
 *  be \c sign * 1.0. 
 * This routine will not catch malformatted numbers like .e10 !
 */
static Real readValue(char*& pos)
{
   assert(isValue(pos));

   char        tmp[MAX_LINE_LEN];
   const char* s = pos;
   char*       t;
   Real        value = 1.0;
   bool        has_digits = false;
   bool        has_emptyexponent = false;

   // 1. sign 
   if ((*s == '+') || (*s == '-'))
      s++;

   // 2. Digits before the decimal dot
   while((*s >= '0') && (*s <= '9'))
   {
      has_digits = true;
      s++;
   }
   // 3. Decimal dot
   if (*s == '.')
   {
      s++;

      // 4. If there was a dot, possible digit behind it
      while((*s >= '0') && (*s <= '9'))
      {
         has_digits = true;
         s++;
      }
   }
   // 5. Exponent
   if (tolower(*s) == 'e')
   {
      has_emptyexponent = true;
      s++;

      // 6. Exponent sign 
      if ((*s == '+') || (*s == '-'))
         s++;

      // 7. Exponent digits
      while((*s >= '0') && (*s <= '9'))
      {
         has_emptyexponent = false;
         s++;      
      }
   }
   assert(s != pos);

   if( has_emptyexponent )
   {
      MSG_WARNING( spxout << "WLPFRD01 Warning: found empty exponent in LP file - check for forbidden variable names with initial 'e' or 'E'\n"; )
   }
   
   if (!has_digits)
      value = (*pos == '-') ? -1.0 : 1.0;
   else
   {
      for(t = tmp; pos != s; pos++)
         *t++ = *pos;   
      *t = '\0';
      value = atof(tmp);
   }
   pos += s - pos;

   assert(pos == s);

   MSG_DEBUG( spxout << "DLPFRD01 readValue = " << value << std::endl; )

   if (isSpace(*pos))
      pos++;

   return value;
}

/// Read the next column name from the input.
/** The name read is looked up and if not found \p emptycol
 *  is added to \p colset. \p pos is advanced behind the name.
 *  @return The Index of the named column.
 */
static int readColName(
   char*& pos, NameSet* colnames, LPColSet& colset, const LPCol* emptycol)
{
   assert(isColName(pos));
   assert(colnames != 0);

   char        name[MAX_LINE_LEN];
   const char* s = pos;
   int         i;
   int         colidx;

   // This are the characters that are not allowed in a column name.
   while((strchr("+-.<>= ", *s) == 0) && (*s != '\0'))
      s++;

   for(i = 0; pos != s; i++, pos++)
      name[i] = *pos;

   name[i] = '\0';

   if ((colidx = colnames->number(name)) < 0)
   {
      // We only add the name if we got an empty column.
      if (emptycol == 0)
         MSG_WARNING( spxout << "WLPFRD02 Unknown variable \"" << name << "\" "; )
      else
      {
         colidx = colnames->num();
         colnames->add(name);
         colset.add(*emptycol);
      }
   }
   MSG_DEBUG( spxout << "DLPFRD03 readColName [" << name << "] = "
                     << colidx << std::endl; )

   if (isSpace(*pos))
      pos++;

   return colidx;
}

/// Read the next <,>,=,==,<=,=<,>=,=> and advance \p pos.
static int readSense(char*& pos)
{
   assert(isSense(pos));

   int sense = *pos++;

   if ((*pos == '<') || (*pos == '>'))
      sense = *pos++;
   else if (*pos == '=')
      pos++;

   MSG_DEBUG( spxout << "DLPFRD04 readSense = " << static_cast<char>(sense)
                     << std::endl; )

   if (isSpace(*pos))
      pos++;

   return sense;
}

/// Is the \p keyword present in \p buf ? If yes, advance \p pos.
/** \p keyword should be lower case. It can contain optional sections
 *  which are enclosed in '[' ']' like "min[imize]". 
 */
static bool hasKeyword(char*& pos, const char* keyword)
{
   int i;
   int k;

   assert(keyword != 0);

   for(i = 0, k = 0; keyword[i] != '\0'; i++, k++)
   {
      if (keyword[i] == '[')
      {
         i++;

         // Here we assumed that we have a ']' for the '['.
         while((tolower(pos[k]) == keyword[i]) && (pos[k] != '\0'))
         {
           k++;
           i++;
         }
         while(keyword[i] != ']')
            i++;         
         --k;
      }
      else
      {
         if (keyword[i] != tolower(pos[k]))
            break;
      }
   }
   // we have to be at the end of the keyword and the word 
   // found on the line also has to end here.
   // Attention: The isSense is a kludge to allow hasKeyword also to process
   //            Inf[inity] keywords in the bounds section.
   if (keyword[i] == '\0' && (pos[k] == '\0' || isSpace(pos[k]) || isSense(&pos[k])))
   {
      pos += k;

      MSG_DEBUG( spxout << "DLPFRD05 hasKeyword: " << keyword << std::endl; )
      return true;
   }
   return false;
}

/// If \p buf start with "name:" store the name in \p rownames 
/// and advance \p pos. 
static bool hasRowName(char*& pos, NameSet* rownames)
{
   const char* s = strchr(pos, ':');

   if (s == 0)
      return false;

   int dcolpos = int(s - pos);

   int end, srt;

   for(end = dcolpos-1; end >= 0; end--)  // skip spaces between name and ":"
     if( pos[end] != ' ')
       break;

   if( end < 0 )  // are there only spaces in front of the ":" ?
   {
     pos = &(pos[dcolpos+1]);
     return false;
   }

   for(srt = end-1; srt >= 0; srt--) // skip spaces in front of name
      if (pos[srt] == ' ')
         break;

   srt++; // go back to the non-space character

   assert( srt <= end && pos[srt] != ' ' );

   char name[MAX_LINE_LEN]; 
   int  i;
   int  k = 0;

   for(i = srt; i <= end; i++)
      name[k++] = pos[i];

   name[k] = '\0';

   if (rownames != 0)
      rownames->add(name);

   pos = &(pos[dcolpos+1]);

   return true;
}

static Real readInfinity(char*& pos)
{
   assert(isInfinity(pos));

   Real sense = (*pos == '-') ? -1.0 : 1.0;

   (void) hasKeyword(++pos, "inf[inity]");

   return sense * infinity;
}

/// Read LP in "CPLEX LP File Format".
/** 
 *  The specification is taken from the
 *  ILOG CPLEX 7.0 Reference Manual, Appendix E, Page 527.
 *
 *  This routine should read (most?) valid LP format files. 
 *  What it will not do, is find all cases where a file is ill formed. 
 *  If this happens it may complain and read nothing or read "something".
 *
 *  Problem: A line ending in '+' or '-' followed by a line starting
 *  with a number, will be regarded as an error.
 * 
 *  The reader will accept the keyword INT[egers] as a synonym for 
 *  GEN[erals] which is an undocumented feature in CPLEX.
 *
 *  A difference to the CPLEX reader, is that no name for the objective 
 *  row is required.
 *
 *  @return true if the file was read correctly
 */  
bool SPxLP::readLPF(
   std::istream& p_input,                ///< input stream.
   NameSet*      p_rnames,               ///< row names.
   NameSet*      p_cnames,               ///< column names.
   DIdxSet*      p_intvars)              ///< integer variables.
{
   enum 
   { 
      START, OBJECTIVE, CONSTRAINTS, BOUNDS, INTEGERS, BINARIES 
   } section = START;

   NameSet*  rnames;                ///< row names.
   NameSet*  cnames;                ///< column names.

   LPColSet  cset;                  ///< the set of columns read.
   LPRowSet  rset;                  ///< the set of rows read.
   LPCol     emptycol;              ///< reusable empty column.
   LPRow     row;                   ///< last assembled row.
   DSVector  vec;                   ///< last assembled vector (from row).
   Real      val = 1.0;
   int       colidx;
   int       sense = 0;

   char      buf [MAX_LINE_LEN];
   char      tmp [MAX_LINE_LEN];
   char      line[MAX_LINE_LEN];
   int       lineno     = 0;
   bool      unnamed    = true;
   bool      finished   = false;
   bool      other;
   bool      have_value = true;
   int       i;
   int       k;
   char*     s;
   char*     pos;
   char*     pos_old = 0;

   cnames = (p_cnames != 0) 
      ? p_cnames : new NameSet();

   cnames->clear();

   try
   {
      rnames = (p_rnames != 0)
         ? p_rnames : new NameSet();
   }catch(std::bad_alloc& x)
   {
      if(p_cnames == 0)
         delete cnames;
      throw x;
   }
   rnames->clear();

   SPxLP::clear(); // clear the LP.

   //--------------------------------------------------------------------------
   //--- Main Loop
   //--------------------------------------------------------------------------
   for(;;)
   {      
      // 0. Read a line from the file.
      if (!p_input.getline(buf, sizeof(buf)))
      {
         if (strlen(buf) == MAX_LINE_LEN - 1)
         {
            MSG_ERROR( spxout << "ELPFRD06 Line exceeds " << MAX_LINE_LEN - 2 
                            << " characters" << std::endl; )
         }
         else
         {
            MSG_ERROR( spxout << "ELPFRD07 No 'End' marker found" << std::endl; )
            finished = true;
         }
         break;
      }
      lineno++;
      i   = 0;
      pos = buf;

      MSG_DEBUG( spxout << "DLPFRD08 Reading line " << lineno
                        << " (pos=" << pos << ")" << std::endl; )

      // 1. Remove comments.
      if (0 != (s = strchr(buf, '\\')))
         *s = '\0';

      // 2. Look for keywords. 
      if (section == START)
      {
         if (hasKeyword(pos, "max[imize]"))
         {
            changeSense(SPxLP::MAXIMIZE);
            section = OBJECTIVE;
         } 
         else if (hasKeyword(pos, "min[imize]"))
         {
            changeSense(SPxLP::MINIMIZE);
            section = OBJECTIVE;
         } 
      }
      else if (section == OBJECTIVE)
      {
         if (hasKeyword(pos, "s[ubject][   ]t[o]")
            || hasKeyword(pos, "s[uch][    ]t[hat]")
            || hasKeyword(pos, "s[.][    ]t[.]")
            || hasKeyword(pos, "lazy con[straints]"))
         {
            // store objective vector            
            for(int j = vec.size() - 1; j >= 0; --j)
               cset.maxObj_w(vec.index(j)) = vec.value(j);
            // multiplication with -1 for minimization is done below
            vec.clear();
            have_value = true;
            val        = 1.0;
            section    = CONSTRAINTS;
         }
      }
      else if (section == CONSTRAINTS &&
              (hasKeyword(pos, "s[ubject][   ]t[o]")
              || hasKeyword(pos, "s[uch][    ]t[hat]")
              || hasKeyword(pos, "s[.][    ]t[.]")))
      {
         have_value = true;
         val        = 1.0;
      }
      else 
      {
         if (hasKeyword(pos, "lazy con[straints]"))
            ;
         else if (hasKeyword(pos, "bound[s]"))
            section = BOUNDS;
         else if (hasKeyword(pos, "bin[ary]"))
            section = BINARIES;
         else if (hasKeyword(pos, "bin[aries]"))
            section = BINARIES;
         else if (hasKeyword(pos, "gen[erals]"))
            section = INTEGERS;
         else if (hasKeyword(pos, "int[egers]")) // this is undocumented
            section = INTEGERS;
         else if (hasKeyword(pos, "end"))
         {
            finished = true;
            break;
         }
         else if (hasKeyword(pos, "s[ubject][   ]t[o]") // second time
            || hasKeyword(pos, "s[uch][    ]t[hat]")
            || hasKeyword(pos, "s[.][    ]t[.]")
            || hasKeyword(pos, "lazy con[straints]"))
         {
            // In principle this has to checked for all keywords above,
            // otherwise we just ignore any half finished constraint
            if (have_value)
               goto syntax_error;

            have_value = true;
            val        = 1.0;
         }
      }

      // 3a. Look for row names in objective and drop it.
      if (section == OBJECTIVE)
         hasRowName(pos, 0);

      // 3b. Look for row name in constraint and store it.
      if (section == CONSTRAINTS)
         if (hasRowName(pos, rnames))
            unnamed = false;

      // 4a. Remove initial spaces.
      while(isSpace(pos[i]))
         i++;

      // 4b. remove spaces if they do not appear before the name of a vaiable.
      for(k = 0; pos[i] != '\0'; i++)
         if (!isSpace(pos[i]) || isColName(&pos[i + 1]))
            tmp[k++] = pos[i];

      tmp[k] = '\0';

      // 5. Is this a empty line ?
      if (tmp[0] == '\0')
         continue;

      // 6. Collapse sequences of '+' and '-'. e.g ++---+ => -
      for(i = 0, k = 0; tmp[i] != '\0'; i++)
      {
         while(((tmp[i] == '+') || (tmp[i] == '-')) 
            && ((tmp[i + 1] == '+') || (tmp[i + 1] == '-')))
         {
            if (tmp[i++] == '-')
               tmp[i] = (tmp[i] == '-') ? '+' : '-';
         }
         line[k++] = tmp[i];
      }
      line[k] = '\0';

      //-----------------------------------------------------------------------
      //--- Line processing loop
      //-----------------------------------------------------------------------
      pos = line;
      
      MSG_DEBUG( spxout << "DLPFRD09 pos=" << pos << std::endl; )

      // 7. We have something left to process. 
      while((pos != 0) && (*pos != '\0'))
      {
         // remember our position, so we are sure we make progress.
         pos_old = pos;

         // Now process the sections 
         switch(section)
         {
         case OBJECTIVE :
            if (isValue(pos))
            {
               Real pre_sign = 1.0;

               /* Allready having here a value could only result from
                * being the first number in a constraint, or a sign
                * '+' or '-' as last token on the previous line.
                */
               if (have_value)
               {
                  if (NE(fabs(val), 1.0))
                     goto syntax_error;
               
                  if (EQ(val, -1.0))
                     pre_sign = val;
               }
               have_value = true;
               val        = readValue(pos) * pre_sign;
            }
            if (*pos == '\0')
               continue;

            if (!have_value || !isColName(pos))
               goto syntax_error;
            
            have_value = false;
            colidx     = readColName(pos, cnames, cset, &emptycol);
            vec.add(colidx, val);
            break;
         case CONSTRAINTS :
            if (isValue(pos))
            {
               Real pre_sign = 1.0;

               /* Allready having here a value could only result from
                * being the first number in a constraint, or a sign
                * '+' or '-' as last token on the previous line.
                */
               if (have_value)
               {
                  if (NE(fabs(val), 1.0))
                     goto syntax_error;
               
                  if (EQ(val, -1.0))
                     pre_sign = val;
               }
                     
               have_value = true;
               val        = readValue(pos) * pre_sign;
               
               if (sense != 0)
               {
                  if (sense == '<')
                  { 
                     row.setLhs(-infinity); 
                     row.setRhs(val);
                  }
                  else if (sense == '>')
                  {
                     row.setLhs(val);
                     row.setRhs(infinity);
                  }
                  else 
                  {
                     assert(sense == '=');
                  
                     row.setLhs(val);
                     row.setRhs(val);
                  }
                  row.setRowVector(vec);
                  rset.add(row);
                  vec.clear();

                  if (!unnamed)
                     unnamed = true;
                  else
                  {
                     char name[16];
                     sprintf(name, "C%d", rset.num());
                     rnames->add(name);
                  }
                  have_value = true;
                  val        = 1.0;
                  sense      = 0;
                  pos        = 0;
                  // next line
                  continue;
               }         
            }
            if (*pos == '\0')
               continue;

            if (have_value)
            {
               if (isColName(pos))
               {
                  colidx = readColName(pos, cnames, cset, &emptycol);

                  if (val != 0.0)
                  {
                     // Do we have this index allready in the row ?
                     int n = vec.number(colidx);

                     // No! So add it.
                     if (n < 0)
                        vec.add(colidx, val);
                     else
                     {                
                        /* Yes. So we add them up and remove the element
                         * if it amounts to zero.
                         */
                        assert(vec.index(n) == colidx);

                        val += vec.value(n);

                        if (val == 0.0)
                           vec.remove(n);
                        else
                           vec.value(n) = val;

                        assert(cnames->has(colidx));

                        MSG_WARNING( spxout << "WLPFRD10 Duplicate index " 
                                            << (*cnames)[colidx] 
                                            << " in line " << lineno 
                                            << std::endl; )
                     }
                  }
                  have_value = false;
               }
               else
               {
                  /* We have a row like c1: <= 5 with no variables in it.
                   * We can not handle 10 <= 5; issue a syntax error.
                   */
                  if (val != 1.0)
                     goto syntax_error;

                  // If the next thing is not the sense we give up also.
                  if (!isSense(pos))
                     goto syntax_error;

                  have_value = false;
               }
            }
            assert(!have_value);

            if (isSense(pos))
               sense = readSense(pos);
            break;
         case BOUNDS :
            other = false;
            sense = 0;

            if (isValue(pos))
            {
               val = isInfinity(pos) ? readInfinity(pos) : readValue(pos);

               if (!isSense(pos))
                  goto syntax_error;

               sense = readSense(pos);
               other = true;
            }
            if (!isColName(pos))
               goto syntax_error;

            if ((colidx = readColName(pos, cnames, cset, 0)) < 0)
            {
               MSG_WARNING( spxout << "WLPFRD11 in Bounds section line " 
                                   << lineno << " ignored" << std::endl; )
                  *pos = '\0';
               continue;
            }
            if (sense)
            {
               if (sense == '<') 
                  cset.lower_w(colidx) = val;
               else if (sense == '>')
                  cset.upper_w(colidx) = val;
               else
               {
                  assert(sense == '=');
                  cset.lower_w(colidx) = val;
                  cset.upper_w(colidx) = val;
               }
            }
            if (isFree(pos))
            {
               cset.lower_w(colidx) = -infinity;
               cset.upper_w(colidx) =  infinity;
               other              = true;
               pos += 4;  // set position after the word "free"
            }
            else if (isSense(pos))
            {
               sense = readSense(pos);
               other = true;

               if (!isValue(pos))
                  goto syntax_error;
               
               val = isInfinity(pos) ? readInfinity(pos) : readValue(pos);

               if (sense == '<') 
                  cset.upper_w(colidx) = val;
               else if (sense == '>')
                  cset.lower_w(colidx) = val;
               else
               {
                  assert(sense == '=');
                  cset.lower_w(colidx) = val;
                  cset.upper_w(colidx) = val;
               }
            }
            /* Do we have only a single column name in the input line?
             * We could ignore this savely, but it is probably a sign 
             * of some other error.
             */
            if (!other)
               goto syntax_error;
            break;
         case BINARIES  :
         case INTEGERS :
            if ((colidx = readColName(pos, cnames, cset, 0)) < 0)
            {
               MSG_WARNING( spxout << "WLPFRD12 in Binary/General section line " 
                                   << lineno << " ignored" << std::endl; )
            }
            else
            {
               if (section == BINARIES)
               {
		  if (cset.lower(colidx) < 0.0)
		  {
		     cset.lower_w(colidx) = 0.0;
		  }
		  if (cset.upper(colidx) > 1.0)
		  {
		     cset.upper_w(colidx) = 1.0;
		  }
               }
               if (p_intvars != 0)
                  p_intvars->addIdx(colidx);
            }
            break;
         case START :
            MSG_ERROR( spxout << "ELPFRD13 This seems to be no LP format file" 
                            << std::endl; )
            goto syntax_error;
         default :
            throw SPxInternalCodeException("XLPFRD01 This should never happen.");
         }
         if (pos == pos_old)
            goto syntax_error;
      }
   }
   //--------------------------------------------------------------------------

   assert(isConsistent());

   addCols(cset);
   assert(isConsistent());

   addRows(rset); 
   assert(isConsistent());

syntax_error:
   if (finished)
   {
      MSG_INFO2( spxout << "ILPFRD14 Finished reading " << lineno
                           << " lines" << std::endl; )
   }
   else
      MSG_ERROR( spxout << "ELPFRD15 Syntax error in line " << lineno << std::endl; )

   if (p_cnames == 0)
      delete cnames;
   if (p_rnames == 0)
      delete rnames;

   MSG_DEBUG( spxout << "DLPFRD16\n" << *this; )

   return finished;
}
} // namespace soplex
