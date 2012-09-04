<?xml version="1.0" encoding="utf-8"?>
<!--
   Copyright (c) 2002 Douglas Gregor <doug.gregor -at- gmail.com>

   Distributed under the Boost Software License, Version 1.0.
   (See accompanying file LICENSE_1_0.txt or copy at
   http://www.boost.org/LICENSE_1_0.txt)
  -->
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
                xmlns:rev="http://www.cs.rpi.edu/~gregod/boost/tools/doc/revision"
                version="1.0">

  <xsl:import href="http://docbook.sourceforge.net/release/xsl/current/htmlhelp/htmlhelp.xsl" />
  <xsl:import href="html-base.xsl" />

  <xsl:param name="generate.toc" />
  <xsl:param name="htmlhelp.enumerate.images" select="1" />
  <xsl:param name="htmlhelp.hhp.tail">
[FILES]
boost.css</xsl:param>

</xsl:stylesheet>
