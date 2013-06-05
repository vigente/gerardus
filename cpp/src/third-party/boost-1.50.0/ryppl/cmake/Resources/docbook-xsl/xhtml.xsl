<?xml version="1.0" encoding="utf-8"?>
<!--
   Copyright (c) 2002 Douglas Gregor <doug.gregor -at- gmail.com>

   Distributed under the Boost Software License, Version 1.0.
   (See accompanying file LICENSE_1_0.txt or copy at
   http://www.boost.org/LICENSE_1_0.txt)
  -->
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
                xmlns="http://www.w3.org/1999/xhtml">

  <!-- Import the HTML chunking stylesheet -->
  <xsl:import href="http://docbook.sourceforge.net/release/xsl/current/xhtml/chunk.xsl" />
  <xsl:import href="http://docbook.sourceforge.net/release/xsl/current/xhtml/math.xsl" />

  <xsl:import href="http://www.boost.org/tools/boostbook/xsl/chunk-common.xsl" />
  <xsl:import href="http://www.boost.org/tools/boostbook/xsl/docbook-layout.xsl" />
  <xsl:import href="http://www.boost.org/tools/boostbook/xsl/admon.xsl" />
  <xsl:import href="http://www.boost.org/tools/boostbook/xsl/xref.xsl" />
  <xsl:import href="http://www.boost.org/tools/boostbook/xsl/relative-href.xsl" />
  <xsl:import href="http://www.boost.org/tools/boostbook/xsl/callout.xsl" />

  <xsl:import href="html-base.xsl" />
  <xsl:import href="navbar.xsl" />

  <xsl:param name="chapter.autolabel" select="1" />
  <xsl:param name="use.id.as.filename" select="1" />

  <xsl:param name="generate.toc">
appendix  toc,title
article/appendix  nop
article   toc,title
book      toc,title
chapter   toc,title
part      toc,title
preface   toc,title
qandadiv  toc
qandaset  toc
reference toc,title
sect1     toc
sect2     toc
sect3     toc
sect4     toc
sect5     toc
section   toc
set       toc,title
  </xsl:param>

  <xsl:template name="output.html.stylesheets">
    <xsl:param name="stylesheets" select="''" />

    <xsl:choose>
      <xsl:when test="contains($stylesheets, ' ')">
        <link rel="stylesheet">
          <xsl:attribute name="href">
            <xsl:call-template name="href.target.relative">
              <xsl:with-param name="target" select="substring-before($stylesheets, ' ')" />
            </xsl:call-template>
          </xsl:attribute>
          <xsl:if test="$html.stylesheet.type != ''">
            <xsl:attribute name="type">
              <xsl:value-of select="$html.stylesheet.type" />
            </xsl:attribute>
          </xsl:if>
        </link>
        <xsl:call-template name="output.html.stylesheets">
          <xsl:with-param name="stylesheets" select="substring-after($stylesheets, ' ')" />
        </xsl:call-template>
      </xsl:when>
      <xsl:when test="$stylesheets != ''">
        <link rel="stylesheet">
          <xsl:attribute name="href">
            <xsl:call-template name="href.target.relative">
              <xsl:with-param name="target" select="$stylesheets" />
            </xsl:call-template>
          </xsl:attribute>
          <xsl:if test="$html.stylesheet.type != ''">
            <xsl:attribute name="type">
              <xsl:value-of select="$html.stylesheet.type" />
            </xsl:attribute>
          </xsl:if>
        </link>
      </xsl:when>
    </xsl:choose>
  </xsl:template>

</xsl:stylesheet>
