<?xml version="1.0" encoding="utf-8"?>
<!--
   Copyright (c) 2002 Douglas Gregor <doug.gregor -at- gmail.com>

   Distributed under the Boost Software License, Version 1.0.
   (See accompanying file LICENSE_1_0.txt or copy at
   http://www.boost.org/LICENSE_1_0.txt)
  -->
<xsl:stylesheet xmlns="http://www.w3.org/1999/xhtml"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
                version="1.0">

  <xsl:template name="header.navigation">
    <xsl:param name="prev" select="/foo" />
    <xsl:param name="next" select="/foo" />
    <xsl:param name="nav.context" />

    <xsl:call-template name="navbar.spirit">
      <xsl:with-param name="prev" select="$prev" />
      <xsl:with-param name="next" select="$next" />
      <xsl:with-param name="nav.context" select="$nav.context" />
    </xsl:call-template>
  </xsl:template>


  <xsl:template name="footer.navigation">
    <xsl:param name="prev" select="/foo" />
    <xsl:param name="next" select="/foo" />
    <xsl:param name="nav.context" />

    <xsl:call-template name="navbar.spirit">
      <xsl:with-param name="prev" select="$prev" />
      <xsl:with-param name="next" select="$next" />
      <xsl:with-param name="nav.context" select="$nav.context" />
    </xsl:call-template>
  </xsl:template>


  <xsl:template name="navbar.spirit">
    <xsl:param name="prev" select="/foo" />
    <xsl:param name="next" select="/foo" />
    <xsl:param name="nav.context" />

    <xsl:variable name="home" select="/*[1]" />
    <xsl:variable name="up" select="parent::*" />

    <div class="navigation">
      <!-- prev -->
      <xsl:if test="count($prev) > 0">
        <a class="prev" accesskey="p">
          <xsl:attribute name="title">
            <xsl:apply-templates select="$prev" mode="object.title.markup" />
          </xsl:attribute>
          <xsl:attribute name="href">
            <xsl:call-template name="href.target">
              <xsl:with-param name="object" select="$prev" />
            </xsl:call-template>
          </xsl:attribute>
          <xsl:text>Prev</xsl:text>
        </a>
      </xsl:if>
      <!-- up -->
      <xsl:if test="count($up) > 0">
        <a class="up" accesskey="u">
          <xsl:attribute name="title">
            <xsl:apply-templates select="$up" mode="object.title.markup" />
          </xsl:attribute>
          <xsl:attribute name="href">
            <xsl:call-template name="href.target">
              <xsl:with-param name="object" select="$up" />
            </xsl:call-template>
          </xsl:attribute>
          <xsl:text>Up</xsl:text>
        </a>
      </xsl:if>
      <!-- home -->
      <xsl:if test="$home != . or $nav.context = 'toc'">
        <a class="home" accesskey="h">
          <xsl:attribute name="title">
            <xsl:apply-templates select="$home" mode="object.title.markup" />
          </xsl:attribute>
          <xsl:attribute name="href">
            <xsl:call-template name="href.target">
              <xsl:with-param name="object" select="$home" />
            </xsl:call-template>
          </xsl:attribute>
          <xsl:text>Home</xsl:text>
        </a>
      </xsl:if>
      <!-- next -->
      <xsl:if test="count($next) > 0">
        <a class="next" accesskey="n">
          <xsl:attribute name="title">
            <xsl:apply-templates select="$next" mode="object.title.markup"/>
          </xsl:attribute>
          <xsl:attribute name="href">
            <xsl:call-template name="href.target">
              <xsl:with-param name="object" select="$next" />
            </xsl:call-template></xsl:attribute>
          <xsl:text>Next</xsl:text>
        </a>
      </xsl:if>
    </div>
  </xsl:template>

</xsl:stylesheet>
