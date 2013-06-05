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

  <xsl:param name="admon.style" />
  <xsl:param name="admon.graphics" select="1" />
  <xsl:param name="boostbook.verbose" select="0" />
  <xsl:param name="css.decoration" select="0" />
  <xsl:param name="html.stylesheet" select="'ryppl.css'" />
  <xsl:param name="navig.graphics" select="1" />
  <xsl:param name="navig.graphics.extension" select="'.png'" />
  <xsl:param name="chapter.autolabel" select="1" />
  <xsl:param name="refentry.generate.name" select="0" />
  <xsl:param name="refentry.generate.title" select="1" />
  <xsl:param name="make.year.ranges" select="1" />
  <xsl:param name="generate.section.toc.level" select="3" />
  <xsl:param name="doc.standalone">false</xsl:param>
  <xsl:param name="chunker.output.indent">yes</xsl:param>
  <xsl:param name="chunker.output.encoding">US-ASCII</xsl:param>
  <xsl:param name="chunk.quietly" select="not(number($boostbook.verbose))"/>
  <xsl:param name="toc.max.depth">2</xsl:param>
  <xsl:param name="callout.graphics.number.limit" select="15" />

  <xsl:template name="format.revision">
    <xsl:param name="text" />

    <xsl:variable name="seperator">
      <xsl:choose>
        <xsl:when test="contains($text, '/')">
          <xsl:value-of select="'/'" />
        </xsl:when>
        <xsl:otherwise>
          <xsl:value-of select="'-'" />
        </xsl:otherwise>
      </xsl:choose>
    </xsl:variable>

    <!-- Remove the "$Date: " -->
    <xsl:variable name="text.noprefix" select="substring-after($text, '$Date: ')" />

    <!-- Grab the year -->
    <xsl:variable name="year" select="substring-before($text.noprefix, $seperator)" />
    <xsl:variable name="text.noyear" select="substring-after($text.noprefix, $seperator)" />

    <!-- Grab the month -->
    <xsl:variable name="month" select="substring-before($text.noyear, $seperator)" />
    <xsl:variable name="text.nomonth" select="substring-after($text.noyear, $seperator)" />

    <!-- Grab the day -->
    <xsl:variable name="day" select="substring-before($text.nomonth, ' ')" />
    <xsl:variable name="text.noday" select="substring-after($text.nomonth, ' ')" />

    <xsl:variable name="month.name">
      <xsl:choose>
        <xsl:when test="$month=1">January</xsl:when>
        <xsl:when test="$month=2">February</xsl:when>
        <xsl:when test="$month=3">March</xsl:when>
        <xsl:when test="$month=4">April</xsl:when>
        <xsl:when test="$month=5">May</xsl:when>
        <xsl:when test="$month=6">June</xsl:when>
        <xsl:when test="$month=7">July</xsl:when>
        <xsl:when test="$month=8">August</xsl:when>
        <xsl:when test="$month=9">September</xsl:when>
        <xsl:when test="$month=10">October</xsl:when>
        <xsl:when test="$month=11">November</xsl:when>
        <xsl:when test="$month=12">December</xsl:when>
      </xsl:choose>
    </xsl:variable>

    <xsl:value-of select="concat($month.name, ' ', $day, ', ', $year)" />
  </xsl:template>

  <!-- Footer Copyright -->
  <xsl:template match="copyright" mode="boost.footer">
    <p>
      <xsl:call-template name="gentext">
        <xsl:with-param name="key" select="'Copyright'" />
      </xsl:call-template>
      <xsl:call-template name="gentext.space" />
      <xsl:call-template name="dingbat">
        <xsl:with-param name="dingbat">copyright</xsl:with-param>
      </xsl:call-template>
      <xsl:call-template name="gentext.space" />
      <xsl:call-template name="copyright.years">
        <xsl:with-param name="years" select="year" />
        <xsl:with-param name="print.ranges" select="$make.year.ranges" />
        <xsl:with-param name="single.year.ranges" select="$make.single.year.ranges" />
      </xsl:call-template>
      <xsl:call-template name="gentext.space" />
      <strong><xsl:apply-templates select="holder" mode="titlepage.mode" /></strong>
    </p>
  </xsl:template>

  <!-- Footer License -->
  <xsl:template match="legalnotice" mode="boost.footer">
    <xsl:apply-templates select="para" mode="titlepage.mode" />
  </xsl:template>

  <xsl:template name="user.footer.content">
    <hr />
    <div class="copyright">
      <xsl:apply-templates select="ancestor-or-self::*/*/copyright" mode="boost.footer" />
      <xsl:apply-templates select="ancestor-or-self::*/*/legalnotice" mode="boost.footer" />

      <xsl:variable name="revision-nodes" select="ancestor-or-self::*[not (attribute::rev:last-revision='')]" />
      <xsl:if test="count($revision-nodes) &gt; 0">
        <xsl:variable name="revision-node" select="$revision-nodes[last()]" />
        <xsl:variable name="revision-text">
          <xsl:value-of select="normalize-space($revision-node/attribute::rev:last-revision)" />
        </xsl:variable>
        <xsl:if test="string-length($revision-text) &gt; 0">
          <p>
            <xsl:text>Last revised: </xsl:text>
            <xsl:call-template name="format.revision">
              <xsl:with-param name="text" select="$revision-text" />
            </xsl:call-template>
          </p>
        </xsl:if>
      </xsl:if>
    </div>
  </xsl:template>

  <!-- We don't want refentry's to show up in the TOC because they will merely be redundant with the synopsis. -->
  <xsl:template match="refentry" mode="toc" />

  <!-- override the behaviour of some DocBook elements for better
       rendering facilities -->

  <xsl:template match = "programlisting[ancestor::informaltable]">
    <pre class = "table-{name(.)}"><xsl:apply-templates /></pre>
  </xsl:template>

  <xsl:template match="refsynopsisdiv">
    <h2 class="{name(.)}-title">Synopsis</h2>
    <div class="{name(.)}">
      <xsl:apply-templates />
    </div>
  </xsl:template>

  <xsl:template name="generate.html.title"/>
  
<!-- ============================================================ -->

  <xsl:template match="author|copyright|legalnotice" mode="titlepage.mode" />

</xsl:stylesheet>
