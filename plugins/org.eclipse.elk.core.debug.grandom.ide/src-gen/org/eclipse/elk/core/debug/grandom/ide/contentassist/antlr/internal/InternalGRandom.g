/*******************************************************************************
 * Copyright (c) 2018 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
grammar InternalGRandom;

options {
	superClass=AbstractInternalContentAssistParser;
}

@lexer::header {
package org.eclipse.elk.core.debug.grandom.ide.contentassist.antlr.internal;

// Hack: Use our own Lexer superclass by means of import. 
// Currently there is no other way to specify the superclass for the lexer.
import org.eclipse.xtext.ide.editor.contentassist.antlr.internal.Lexer;
}

@parser::header {
package org.eclipse.elk.core.debug.grandom.ide.contentassist.antlr.internal;

import java.io.InputStream;
import org.eclipse.xtext.*;
import org.eclipse.xtext.parser.*;
import org.eclipse.xtext.parser.impl.*;
import org.eclipse.emf.ecore.util.EcoreUtil;
import org.eclipse.emf.ecore.EObject;
import org.eclipse.xtext.parser.antlr.XtextTokenStream;
import org.eclipse.xtext.parser.antlr.XtextTokenStream.HiddenTokens;
import org.eclipse.xtext.ide.editor.contentassist.antlr.internal.AbstractInternalContentAssistParser;
import org.eclipse.xtext.ide.editor.contentassist.antlr.internal.DFA;
import org.eclipse.elk.core.debug.grandom.services.GRandomGrammarAccess;

}
@parser::members {
	private GRandomGrammarAccess grammarAccess;

	public void setGrammarAccess(GRandomGrammarAccess grammarAccess) {
		this.grammarAccess = grammarAccess;
	}

	@Override
	protected Grammar getGrammar() {
		return grammarAccess.getGrammar();
	}

	@Override
	protected String getValueForTokenName(String tokenName) {
		return tokenName;
	}
}

// Entry rule entryRuleRandGraph
entryRuleRandGraph
:
{ before(grammarAccess.getRandGraphRule()); }
	 ruleRandGraph
{ after(grammarAccess.getRandGraphRule()); } 
	 EOF 
;

// Rule RandGraph
ruleRandGraph 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getRandGraphAccess().getConfigsAssignment()); }
		(rule__RandGraph__ConfigsAssignment)*
		{ after(grammarAccess.getRandGraphAccess().getConfigsAssignment()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

// Entry rule entryRuleConfiguration
entryRuleConfiguration
:
{ before(grammarAccess.getConfigurationRule()); }
	 ruleConfiguration
{ after(grammarAccess.getConfigurationRule()); } 
	 EOF 
;

// Rule Configuration
ruleConfiguration 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getConfigurationAccess().getGroup()); }
		(rule__Configuration__Group__0)
		{ after(grammarAccess.getConfigurationAccess().getGroup()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

// Entry rule entryRuleHierarchy
entryRuleHierarchy
:
{ before(grammarAccess.getHierarchyRule()); }
	 ruleHierarchy
{ after(grammarAccess.getHierarchyRule()); } 
	 EOF 
;

// Rule Hierarchy
ruleHierarchy 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getHierarchyAccess().getGroup()); }
		(rule__Hierarchy__Group__0)
		{ after(grammarAccess.getHierarchyAccess().getGroup()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

// Entry rule entryRuleEdges
entryRuleEdges
:
{ before(grammarAccess.getEdgesRule()); }
	 ruleEdges
{ after(grammarAccess.getEdgesRule()); } 
	 EOF 
;

// Rule Edges
ruleEdges 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getEdgesAccess().getGroup()); }
		(rule__Edges__Group__0)
		{ after(grammarAccess.getEdgesAccess().getGroup()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

// Entry rule entryRuleNodes
entryRuleNodes
:
{ before(grammarAccess.getNodesRule()); }
	 ruleNodes
{ after(grammarAccess.getNodesRule()); } 
	 EOF 
;

// Rule Nodes
ruleNodes 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getNodesAccess().getGroup()); }
		(rule__Nodes__Group__0)
		{ after(grammarAccess.getNodesAccess().getGroup()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

// Entry rule entryRuleSize
entryRuleSize
:
{ before(grammarAccess.getSizeRule()); }
	 ruleSize
{ after(grammarAccess.getSizeRule()); } 
	 EOF 
;

// Rule Size
ruleSize 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getSizeAccess().getGroup()); }
		(rule__Size__Group__0)
		{ after(grammarAccess.getSizeAccess().getGroup()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

// Entry rule entryRulePorts
entryRulePorts
:
{ before(grammarAccess.getPortsRule()); }
	 rulePorts
{ after(grammarAccess.getPortsRule()); } 
	 EOF 
;

// Rule Ports
rulePorts 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getPortsAccess().getGroup()); }
		(rule__Ports__Group__0)
		{ after(grammarAccess.getPortsAccess().getGroup()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

// Entry rule entryRuleFlow
entryRuleFlow
:
{ before(grammarAccess.getFlowRule()); }
	 ruleFlow
{ after(grammarAccess.getFlowRule()); } 
	 EOF 
;

// Rule Flow
ruleFlow 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getFlowAccess().getGroup()); }
		(rule__Flow__Group__0)
		{ after(grammarAccess.getFlowAccess().getGroup()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

// Entry rule entryRuleLabels
entryRuleLabels
:
{ before(grammarAccess.getLabelsRule()); }
	 ruleLabels
{ after(grammarAccess.getLabelsRule()); } 
	 EOF 
;

// Rule Labels
ruleLabels 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getLabelsAccess().getLabelsKeyword()); }
		'labels'
		{ after(grammarAccess.getLabelsAccess().getLabelsKeyword()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

// Entry rule entryRuleDoubleQuantity
entryRuleDoubleQuantity
:
{ before(grammarAccess.getDoubleQuantityRule()); }
	 ruleDoubleQuantity
{ after(grammarAccess.getDoubleQuantityRule()); } 
	 EOF 
;

// Rule DoubleQuantity
ruleDoubleQuantity 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getDoubleQuantityAccess().getAlternatives()); }
		(rule__DoubleQuantity__Alternatives)
		{ after(grammarAccess.getDoubleQuantityAccess().getAlternatives()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

// Entry rule entryRulePm
entryRulePm
:
{ before(grammarAccess.getPmRule()); }
	 rulePm
{ after(grammarAccess.getPmRule()); } 
	 EOF 
;

// Rule Pm
rulePm 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getPmAccess().getPlusSignSolidusHyphenMinusKeyword()); }
		'+/-'
		{ after(grammarAccess.getPmAccess().getPlusSignSolidusHyphenMinusKeyword()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

// Entry rule entryRuleDouble
entryRuleDouble
:
{ before(grammarAccess.getDoubleRule()); }
	 ruleDouble
{ after(grammarAccess.getDoubleRule()); } 
	 EOF 
;

// Rule Double
ruleDouble 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getDoubleAccess().getGroup()); }
		(rule__Double__Group__0)
		{ after(grammarAccess.getDoubleAccess().getGroup()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

// Entry rule entryRuleInteger
entryRuleInteger
:
{ before(grammarAccess.getIntegerRule()); }
	 ruleInteger
{ after(grammarAccess.getIntegerRule()); } 
	 EOF 
;

// Rule Integer
ruleInteger 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getIntegerAccess().getGroup()); }
		(rule__Integer__Group__0)
		{ after(grammarAccess.getIntegerAccess().getGroup()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

// Rule Formats
ruleFormats
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getFormatsAccess().getAlternatives()); }
		(rule__Formats__Alternatives)
		{ after(grammarAccess.getFormatsAccess().getAlternatives()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

// Rule Form
ruleForm
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getFormAccess().getAlternatives()); }
		(rule__Form__Alternatives)
		{ after(grammarAccess.getFormAccess().getAlternatives()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

// Rule Side
ruleSide
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getSideAccess().getAlternatives()); }
		(rule__Side__Alternatives)
		{ after(grammarAccess.getSideAccess().getAlternatives()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

// Rule FlowType
ruleFlowType
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getFlowTypeAccess().getAlternatives()); }
		(rule__FlowType__Alternatives)
		{ after(grammarAccess.getFlowTypeAccess().getAlternatives()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

// Rule ConstraintType
ruleConstraintType
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getConstraintTypeAccess().getAlternatives()); }
		(rule__ConstraintType__Alternatives)
		{ after(grammarAccess.getConstraintTypeAccess().getAlternatives()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__Edges__Alternatives_0_1
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getEdgesAccess().getDensityAssignment_0_1_0()); }
		(rule__Edges__DensityAssignment_0_1_0)
		{ after(grammarAccess.getEdgesAccess().getDensityAssignment_0_1_0()); }
	)
	|
	(
		{ before(grammarAccess.getEdgesAccess().getTotalAssignment_0_1_1()); }
		(rule__Edges__TotalAssignment_0_1_1)
		{ after(grammarAccess.getEdgesAccess().getTotalAssignment_0_1_1()); }
	)
	|
	(
		{ before(grammarAccess.getEdgesAccess().getRelativeAssignment_0_1_2()); }
		(rule__Edges__RelativeAssignment_0_1_2)
		{ after(grammarAccess.getEdgesAccess().getRelativeAssignment_0_1_2()); }
	)
	|
	(
		{ before(grammarAccess.getEdgesAccess().getOutboundAssignment_0_1_3()); }
		(rule__Edges__OutboundAssignment_0_1_3)
		{ after(grammarAccess.getEdgesAccess().getOutboundAssignment_0_1_3()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__DoubleQuantity__Alternatives
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getDoubleQuantityAccess().getQuantAssignment_0()); }
		(rule__DoubleQuantity__QuantAssignment_0)
		{ after(grammarAccess.getDoubleQuantityAccess().getQuantAssignment_0()); }
	)
	|
	(
		{ before(grammarAccess.getDoubleQuantityAccess().getGroup_1()); }
		(rule__DoubleQuantity__Group_1__0)
		{ after(grammarAccess.getDoubleQuantityAccess().getGroup_1()); }
	)
	|
	(
		{ before(grammarAccess.getDoubleQuantityAccess().getGroup_2()); }
		(rule__DoubleQuantity__Group_2__0)
		{ after(grammarAccess.getDoubleQuantityAccess().getGroup_2()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__Formats__Alternatives
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getFormatsAccess().getElktEnumLiteralDeclaration_0()); }
		('elkt')
		{ after(grammarAccess.getFormatsAccess().getElktEnumLiteralDeclaration_0()); }
	)
	|
	(
		{ before(grammarAccess.getFormatsAccess().getElkgEnumLiteralDeclaration_1()); }
		('elkg')
		{ after(grammarAccess.getFormatsAccess().getElkgEnumLiteralDeclaration_1()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__Form__Alternatives
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getFormAccess().getTreesEnumLiteralDeclaration_0()); }
		('trees')
		{ after(grammarAccess.getFormAccess().getTreesEnumLiteralDeclaration_0()); }
	)
	|
	(
		{ before(grammarAccess.getFormAccess().getCustomEnumLiteralDeclaration_1()); }
		('graphs')
		{ after(grammarAccess.getFormAccess().getCustomEnumLiteralDeclaration_1()); }
	)
	|
	(
		{ before(grammarAccess.getFormAccess().getBipartiteEnumLiteralDeclaration_2()); }
		('bipartite graphs')
		{ after(grammarAccess.getFormAccess().getBipartiteEnumLiteralDeclaration_2()); }
	)
	|
	(
		{ before(grammarAccess.getFormAccess().getBiconnectedEnumLiteralDeclaration_3()); }
		('biconnected graphs')
		{ after(grammarAccess.getFormAccess().getBiconnectedEnumLiteralDeclaration_3()); }
	)
	|
	(
		{ before(grammarAccess.getFormAccess().getTriconnectedEnumLiteralDeclaration_4()); }
		('triconnected graphs')
		{ after(grammarAccess.getFormAccess().getTriconnectedEnumLiteralDeclaration_4()); }
	)
	|
	(
		{ before(grammarAccess.getFormAccess().getAcyclicEnumLiteralDeclaration_5()); }
		('acyclic graphs')
		{ after(grammarAccess.getFormAccess().getAcyclicEnumLiteralDeclaration_5()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__Side__Alternatives
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getSideAccess().getNorthEnumLiteralDeclaration_0()); }
		('north')
		{ after(grammarAccess.getSideAccess().getNorthEnumLiteralDeclaration_0()); }
	)
	|
	(
		{ before(grammarAccess.getSideAccess().getEastEnumLiteralDeclaration_1()); }
		('east')
		{ after(grammarAccess.getSideAccess().getEastEnumLiteralDeclaration_1()); }
	)
	|
	(
		{ before(grammarAccess.getSideAccess().getSouthEnumLiteralDeclaration_2()); }
		('south')
		{ after(grammarAccess.getSideAccess().getSouthEnumLiteralDeclaration_2()); }
	)
	|
	(
		{ before(grammarAccess.getSideAccess().getWestEnumLiteralDeclaration_3()); }
		('west')
		{ after(grammarAccess.getSideAccess().getWestEnumLiteralDeclaration_3()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__FlowType__Alternatives
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getFlowTypeAccess().getIncomingEnumLiteralDeclaration_0()); }
		('incoming')
		{ after(grammarAccess.getFlowTypeAccess().getIncomingEnumLiteralDeclaration_0()); }
	)
	|
	(
		{ before(grammarAccess.getFlowTypeAccess().getOutgoingEnumLiteralDeclaration_1()); }
		('outgoing')
		{ after(grammarAccess.getFlowTypeAccess().getOutgoingEnumLiteralDeclaration_1()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ConstraintType__Alternatives
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getConstraintTypeAccess().getFreeEnumLiteralDeclaration_0()); }
		('free')
		{ after(grammarAccess.getConstraintTypeAccess().getFreeEnumLiteralDeclaration_0()); }
	)
	|
	(
		{ before(grammarAccess.getConstraintTypeAccess().getSideEnumLiteralDeclaration_1()); }
		('side')
		{ after(grammarAccess.getConstraintTypeAccess().getSideEnumLiteralDeclaration_1()); }
	)
	|
	(
		{ before(grammarAccess.getConstraintTypeAccess().getPositionEnumLiteralDeclaration_2()); }
		('position')
		{ after(grammarAccess.getConstraintTypeAccess().getPositionEnumLiteralDeclaration_2()); }
	)
	|
	(
		{ before(grammarAccess.getConstraintTypeAccess().getOrderEnumLiteralDeclaration_3()); }
		('order')
		{ after(grammarAccess.getConstraintTypeAccess().getOrderEnumLiteralDeclaration_3()); }
	)
	|
	(
		{ before(grammarAccess.getConstraintTypeAccess().getRatioEnumLiteralDeclaration_4()); }
		('ratio')
		{ after(grammarAccess.getConstraintTypeAccess().getRatioEnumLiteralDeclaration_4()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__Group__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Configuration__Group__0__Impl
	rule__Configuration__Group__1
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__Group__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getConfigurationAccess().getGenerateKeyword_0()); }
	'generate'
	{ after(grammarAccess.getConfigurationAccess().getGenerateKeyword_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__Group__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Configuration__Group__1__Impl
	rule__Configuration__Group__2
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__Group__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getConfigurationAccess().getSamplesAssignment_1()); }
	(rule__Configuration__SamplesAssignment_1)
	{ after(grammarAccess.getConfigurationAccess().getSamplesAssignment_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__Group__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Configuration__Group__2__Impl
	rule__Configuration__Group__3
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__Group__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getConfigurationAccess().getFormAssignment_2()); }
	(rule__Configuration__FormAssignment_2)
	{ after(grammarAccess.getConfigurationAccess().getFormAssignment_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__Group__3
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Configuration__Group__3__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__Group__3__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getConfigurationAccess().getGroup_3()); }
	(rule__Configuration__Group_3__0)?
	{ after(grammarAccess.getConfigurationAccess().getGroup_3()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__Configuration__Group_3__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Configuration__Group_3__0__Impl
	rule__Configuration__Group_3__1
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__Group_3__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getConfigurationAccess().getLeftCurlyBracketKeyword_3_0()); }
	'{'
	{ after(grammarAccess.getConfigurationAccess().getLeftCurlyBracketKeyword_3_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__Group_3__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Configuration__Group_3__1__Impl
	rule__Configuration__Group_3__2
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__Group_3__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getConfigurationAccess().getUnorderedGroup_3_1()); }
	(rule__Configuration__UnorderedGroup_3_1)
	{ after(grammarAccess.getConfigurationAccess().getUnorderedGroup_3_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__Group_3__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Configuration__Group_3__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__Group_3__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getConfigurationAccess().getRightCurlyBracketKeyword_3_2()); }
	'}'
	{ after(grammarAccess.getConfigurationAccess().getRightCurlyBracketKeyword_3_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__Configuration__Group_3_1_2__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Configuration__Group_3_1_2__0__Impl
	rule__Configuration__Group_3_1_2__1
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__Group_3_1_2__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getConfigurationAccess().getMWAssignment_3_1_2_0()); }
	(rule__Configuration__MWAssignment_3_1_2_0)
	{ after(grammarAccess.getConfigurationAccess().getMWAssignment_3_1_2_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__Group_3_1_2__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Configuration__Group_3_1_2__1__Impl
	rule__Configuration__Group_3_1_2__2
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__Group_3_1_2__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getConfigurationAccess().getEqualsSignKeyword_3_1_2_1()); }
	'='
	{ after(grammarAccess.getConfigurationAccess().getEqualsSignKeyword_3_1_2_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__Group_3_1_2__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Configuration__Group_3_1_2__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__Group_3_1_2__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getConfigurationAccess().getMaxWidthAssignment_3_1_2_2()); }
	(rule__Configuration__MaxWidthAssignment_3_1_2_2)
	{ after(grammarAccess.getConfigurationAccess().getMaxWidthAssignment_3_1_2_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__Configuration__Group_3_1_3__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Configuration__Group_3_1_3__0__Impl
	rule__Configuration__Group_3_1_3__1
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__Group_3_1_3__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getConfigurationAccess().getMDAssignment_3_1_3_0()); }
	(rule__Configuration__MDAssignment_3_1_3_0)
	{ after(grammarAccess.getConfigurationAccess().getMDAssignment_3_1_3_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__Group_3_1_3__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Configuration__Group_3_1_3__1__Impl
	rule__Configuration__Group_3_1_3__2
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__Group_3_1_3__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getConfigurationAccess().getEqualsSignKeyword_3_1_3_1()); }
	'='
	{ after(grammarAccess.getConfigurationAccess().getEqualsSignKeyword_3_1_3_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__Group_3_1_3__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Configuration__Group_3_1_3__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__Group_3_1_3__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getConfigurationAccess().getMaxDegreeAssignment_3_1_3_2()); }
	(rule__Configuration__MaxDegreeAssignment_3_1_3_2)
	{ after(grammarAccess.getConfigurationAccess().getMaxDegreeAssignment_3_1_3_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__Configuration__Group_3_1_4__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Configuration__Group_3_1_4__0__Impl
	rule__Configuration__Group_3_1_4__1
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__Group_3_1_4__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getConfigurationAccess().getPFAssignment_3_1_4_0()); }
	(rule__Configuration__PFAssignment_3_1_4_0)
	{ after(grammarAccess.getConfigurationAccess().getPFAssignment_3_1_4_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__Group_3_1_4__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Configuration__Group_3_1_4__1__Impl
	rule__Configuration__Group_3_1_4__2
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__Group_3_1_4__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getConfigurationAccess().getEqualsSignKeyword_3_1_4_1()); }
	'='
	{ after(grammarAccess.getConfigurationAccess().getEqualsSignKeyword_3_1_4_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__Group_3_1_4__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Configuration__Group_3_1_4__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__Group_3_1_4__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getConfigurationAccess().getFractionAssignment_3_1_4_2()); }
	(rule__Configuration__FractionAssignment_3_1_4_2)
	{ after(grammarAccess.getConfigurationAccess().getFractionAssignment_3_1_4_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__Configuration__Group_3_1_6__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Configuration__Group_3_1_6__0__Impl
	rule__Configuration__Group_3_1_6__1
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__Group_3_1_6__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getConfigurationAccess().getSeedKeyword_3_1_6_0()); }
	'seed'
	{ after(grammarAccess.getConfigurationAccess().getSeedKeyword_3_1_6_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__Group_3_1_6__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Configuration__Group_3_1_6__1__Impl
	rule__Configuration__Group_3_1_6__2
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__Group_3_1_6__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getConfigurationAccess().getEqualsSignKeyword_3_1_6_1()); }
	'='
	{ after(grammarAccess.getConfigurationAccess().getEqualsSignKeyword_3_1_6_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__Group_3_1_6__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Configuration__Group_3_1_6__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__Group_3_1_6__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getConfigurationAccess().getSeedAssignment_3_1_6_2()); }
	(rule__Configuration__SeedAssignment_3_1_6_2)
	{ after(grammarAccess.getConfigurationAccess().getSeedAssignment_3_1_6_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__Configuration__Group_3_1_7__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Configuration__Group_3_1_7__0__Impl
	rule__Configuration__Group_3_1_7__1
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__Group_3_1_7__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getConfigurationAccess().getFormatKeyword_3_1_7_0()); }
	'format'
	{ after(grammarAccess.getConfigurationAccess().getFormatKeyword_3_1_7_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__Group_3_1_7__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Configuration__Group_3_1_7__1__Impl
	rule__Configuration__Group_3_1_7__2
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__Group_3_1_7__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getConfigurationAccess().getEqualsSignKeyword_3_1_7_1()); }
	'='
	{ after(grammarAccess.getConfigurationAccess().getEqualsSignKeyword_3_1_7_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__Group_3_1_7__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Configuration__Group_3_1_7__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__Group_3_1_7__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getConfigurationAccess().getFormatAssignment_3_1_7_2()); }
	(rule__Configuration__FormatAssignment_3_1_7_2)
	{ after(grammarAccess.getConfigurationAccess().getFormatAssignment_3_1_7_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__Configuration__Group_3_1_8__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Configuration__Group_3_1_8__0__Impl
	rule__Configuration__Group_3_1_8__1
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__Group_3_1_8__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getConfigurationAccess().getFilenameKeyword_3_1_8_0()); }
	'filename'
	{ after(grammarAccess.getConfigurationAccess().getFilenameKeyword_3_1_8_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__Group_3_1_8__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Configuration__Group_3_1_8__1__Impl
	rule__Configuration__Group_3_1_8__2
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__Group_3_1_8__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getConfigurationAccess().getEqualsSignKeyword_3_1_8_1()); }
	'='
	{ after(grammarAccess.getConfigurationAccess().getEqualsSignKeyword_3_1_8_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__Group_3_1_8__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Configuration__Group_3_1_8__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__Group_3_1_8__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getConfigurationAccess().getFilenameAssignment_3_1_8_2()); }
	(rule__Configuration__FilenameAssignment_3_1_8_2)
	{ after(grammarAccess.getConfigurationAccess().getFilenameAssignment_3_1_8_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__Hierarchy__Group__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Hierarchy__Group__0__Impl
	rule__Hierarchy__Group__1
;
finally {
	restoreStackSize(stackSize);
}

rule__Hierarchy__Group__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getHierarchyAccess().getHierarchyAction_0()); }
	()
	{ after(grammarAccess.getHierarchyAccess().getHierarchyAction_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Hierarchy__Group__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Hierarchy__Group__1__Impl
	rule__Hierarchy__Group__2
;
finally {
	restoreStackSize(stackSize);
}

rule__Hierarchy__Group__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getHierarchyAccess().getHierarchyKeyword_1()); }
	'hierarchy'
	{ after(grammarAccess.getHierarchyAccess().getHierarchyKeyword_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Hierarchy__Group__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Hierarchy__Group__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__Hierarchy__Group__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getHierarchyAccess().getGroup_2()); }
	(rule__Hierarchy__Group_2__0)?
	{ after(grammarAccess.getHierarchyAccess().getGroup_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__Hierarchy__Group_2__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Hierarchy__Group_2__0__Impl
	rule__Hierarchy__Group_2__1
;
finally {
	restoreStackSize(stackSize);
}

rule__Hierarchy__Group_2__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getHierarchyAccess().getLeftCurlyBracketKeyword_2_0()); }
	'{'
	{ after(grammarAccess.getHierarchyAccess().getLeftCurlyBracketKeyword_2_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Hierarchy__Group_2__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Hierarchy__Group_2__1__Impl
	rule__Hierarchy__Group_2__2
;
finally {
	restoreStackSize(stackSize);
}

rule__Hierarchy__Group_2__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getHierarchyAccess().getUnorderedGroup_2_1()); }
	(rule__Hierarchy__UnorderedGroup_2_1)
	{ after(grammarAccess.getHierarchyAccess().getUnorderedGroup_2_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Hierarchy__Group_2__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Hierarchy__Group_2__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__Hierarchy__Group_2__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getHierarchyAccess().getRightCurlyBracketKeyword_2_2()); }
	'}'
	{ after(grammarAccess.getHierarchyAccess().getRightCurlyBracketKeyword_2_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__Hierarchy__Group_2_1_0__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Hierarchy__Group_2_1_0__0__Impl
	rule__Hierarchy__Group_2_1_0__1
;
finally {
	restoreStackSize(stackSize);
}

rule__Hierarchy__Group_2_1_0__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getHierarchyAccess().getLevelsKeyword_2_1_0_0()); }
	'levels'
	{ after(grammarAccess.getHierarchyAccess().getLevelsKeyword_2_1_0_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Hierarchy__Group_2_1_0__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Hierarchy__Group_2_1_0__1__Impl
	rule__Hierarchy__Group_2_1_0__2
;
finally {
	restoreStackSize(stackSize);
}

rule__Hierarchy__Group_2_1_0__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getHierarchyAccess().getEqualsSignKeyword_2_1_0_1()); }
	'='
	{ after(grammarAccess.getHierarchyAccess().getEqualsSignKeyword_2_1_0_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Hierarchy__Group_2_1_0__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Hierarchy__Group_2_1_0__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__Hierarchy__Group_2_1_0__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getHierarchyAccess().getLevelsAssignment_2_1_0_2()); }
	(rule__Hierarchy__LevelsAssignment_2_1_0_2)
	{ after(grammarAccess.getHierarchyAccess().getLevelsAssignment_2_1_0_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__Hierarchy__Group_2_1_1__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Hierarchy__Group_2_1_1__0__Impl
	rule__Hierarchy__Group_2_1_1__1
;
finally {
	restoreStackSize(stackSize);
}

rule__Hierarchy__Group_2_1_1__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getHierarchyAccess().getCrossHierarchyEdgesKeyword_2_1_1_0()); }
	'cross-hierarchy edges'
	{ after(grammarAccess.getHierarchyAccess().getCrossHierarchyEdgesKeyword_2_1_1_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Hierarchy__Group_2_1_1__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Hierarchy__Group_2_1_1__1__Impl
	rule__Hierarchy__Group_2_1_1__2
;
finally {
	restoreStackSize(stackSize);
}

rule__Hierarchy__Group_2_1_1__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getHierarchyAccess().getEqualsSignKeyword_2_1_1_1()); }
	'='
	{ after(grammarAccess.getHierarchyAccess().getEqualsSignKeyword_2_1_1_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Hierarchy__Group_2_1_1__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Hierarchy__Group_2_1_1__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__Hierarchy__Group_2_1_1__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getHierarchyAccess().getEdgesAssignment_2_1_1_2()); }
	(rule__Hierarchy__EdgesAssignment_2_1_1_2)
	{ after(grammarAccess.getHierarchyAccess().getEdgesAssignment_2_1_1_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__Hierarchy__Group_2_1_2__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Hierarchy__Group_2_1_2__0__Impl
	rule__Hierarchy__Group_2_1_2__1
;
finally {
	restoreStackSize(stackSize);
}

rule__Hierarchy__Group_2_1_2__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getHierarchyAccess().getCompoundNodesKeyword_2_1_2_0()); }
	'compound nodes'
	{ after(grammarAccess.getHierarchyAccess().getCompoundNodesKeyword_2_1_2_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Hierarchy__Group_2_1_2__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Hierarchy__Group_2_1_2__1__Impl
	rule__Hierarchy__Group_2_1_2__2
;
finally {
	restoreStackSize(stackSize);
}

rule__Hierarchy__Group_2_1_2__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getHierarchyAccess().getEqualsSignKeyword_2_1_2_1()); }
	'='
	{ after(grammarAccess.getHierarchyAccess().getEqualsSignKeyword_2_1_2_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Hierarchy__Group_2_1_2__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Hierarchy__Group_2_1_2__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__Hierarchy__Group_2_1_2__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getHierarchyAccess().getNumHierarchNodesAssignment_2_1_2_2()); }
	(rule__Hierarchy__NumHierarchNodesAssignment_2_1_2_2)
	{ after(grammarAccess.getHierarchyAccess().getNumHierarchNodesAssignment_2_1_2_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__Hierarchy__Group_2_1_3__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Hierarchy__Group_2_1_3__0__Impl
	rule__Hierarchy__Group_2_1_3__1
;
finally {
	restoreStackSize(stackSize);
}

rule__Hierarchy__Group_2_1_3__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getHierarchyAccess().getCrossHierarchyRelativeEdgesKeyword_2_1_3_0()); }
	'cross-hierarchy relative edges'
	{ after(grammarAccess.getHierarchyAccess().getCrossHierarchyRelativeEdgesKeyword_2_1_3_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Hierarchy__Group_2_1_3__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Hierarchy__Group_2_1_3__1__Impl
	rule__Hierarchy__Group_2_1_3__2
;
finally {
	restoreStackSize(stackSize);
}

rule__Hierarchy__Group_2_1_3__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getHierarchyAccess().getEqualsSignKeyword_2_1_3_1()); }
	'='
	{ after(grammarAccess.getHierarchyAccess().getEqualsSignKeyword_2_1_3_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Hierarchy__Group_2_1_3__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Hierarchy__Group_2_1_3__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__Hierarchy__Group_2_1_3__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getHierarchyAccess().getCrossHierarchRelAssignment_2_1_3_2()); }
	(rule__Hierarchy__CrossHierarchRelAssignment_2_1_3_2)
	{ after(grammarAccess.getHierarchyAccess().getCrossHierarchRelAssignment_2_1_3_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__Edges__Group__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Edges__Group__0__Impl
	rule__Edges__Group__1
;
finally {
	restoreStackSize(stackSize);
}

rule__Edges__Group__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getEdgesAccess().getGroup_0()); }
	(rule__Edges__Group_0__0)
	{ after(grammarAccess.getEdgesAccess().getGroup_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Edges__Group__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Edges__Group__1__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__Edges__Group__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getEdgesAccess().getGroup_1()); }
	(rule__Edges__Group_1__0)?
	{ after(grammarAccess.getEdgesAccess().getGroup_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__Edges__Group_0__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Edges__Group_0__0__Impl
	rule__Edges__Group_0__1
;
finally {
	restoreStackSize(stackSize);
}

rule__Edges__Group_0__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getEdgesAccess().getEdgesKeyword_0_0()); }
	'edges'
	{ after(grammarAccess.getEdgesAccess().getEdgesKeyword_0_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Edges__Group_0__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Edges__Group_0__1__Impl
	rule__Edges__Group_0__2
;
finally {
	restoreStackSize(stackSize);
}

rule__Edges__Group_0__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getEdgesAccess().getAlternatives_0_1()); }
	(rule__Edges__Alternatives_0_1)
	{ after(grammarAccess.getEdgesAccess().getAlternatives_0_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Edges__Group_0__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Edges__Group_0__2__Impl
	rule__Edges__Group_0__3
;
finally {
	restoreStackSize(stackSize);
}

rule__Edges__Group_0__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getEdgesAccess().getEqualsSignKeyword_0_2()); }
	'='
	{ after(grammarAccess.getEdgesAccess().getEqualsSignKeyword_0_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Edges__Group_0__3
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Edges__Group_0__3__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__Edges__Group_0__3__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getEdgesAccess().getNEdgesAssignment_0_3()); }
	(rule__Edges__NEdgesAssignment_0_3)
	{ after(grammarAccess.getEdgesAccess().getNEdgesAssignment_0_3()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__Edges__Group_1__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Edges__Group_1__0__Impl
	rule__Edges__Group_1__1
;
finally {
	restoreStackSize(stackSize);
}

rule__Edges__Group_1__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getEdgesAccess().getLeftCurlyBracketKeyword_1_0()); }
	'{'
	{ after(grammarAccess.getEdgesAccess().getLeftCurlyBracketKeyword_1_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Edges__Group_1__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Edges__Group_1__1__Impl
	rule__Edges__Group_1__2
;
finally {
	restoreStackSize(stackSize);
}

rule__Edges__Group_1__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getEdgesAccess().getUnorderedGroup_1_1()); }
	(rule__Edges__UnorderedGroup_1_1)
	{ after(grammarAccess.getEdgesAccess().getUnorderedGroup_1_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Edges__Group_1__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Edges__Group_1__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__Edges__Group_1__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getEdgesAccess().getRightCurlyBracketKeyword_1_2()); }
	'}'
	{ after(grammarAccess.getEdgesAccess().getRightCurlyBracketKeyword_1_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__Nodes__Group__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Nodes__Group__0__Impl
	rule__Nodes__Group__1
;
finally {
	restoreStackSize(stackSize);
}

rule__Nodes__Group__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getNodesAccess().getNodesAction_0()); }
	()
	{ after(grammarAccess.getNodesAccess().getNodesAction_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Nodes__Group__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Nodes__Group__1__Impl
	rule__Nodes__Group__2
;
finally {
	restoreStackSize(stackSize);
}

rule__Nodes__Group__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getNodesAccess().getNodesKeyword_1()); }
	'nodes'
	{ after(grammarAccess.getNodesAccess().getNodesKeyword_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Nodes__Group__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Nodes__Group__2__Impl
	rule__Nodes__Group__3
;
finally {
	restoreStackSize(stackSize);
}

rule__Nodes__Group__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getNodesAccess().getEqualsSignKeyword_2()); }
	'='
	{ after(grammarAccess.getNodesAccess().getEqualsSignKeyword_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Nodes__Group__3
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Nodes__Group__3__Impl
	rule__Nodes__Group__4
;
finally {
	restoreStackSize(stackSize);
}

rule__Nodes__Group__3__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getNodesAccess().getNNodesAssignment_3()); }
	(rule__Nodes__NNodesAssignment_3)
	{ after(grammarAccess.getNodesAccess().getNNodesAssignment_3()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Nodes__Group__4
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Nodes__Group__4__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__Nodes__Group__4__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getNodesAccess().getGroup_4()); }
	(rule__Nodes__Group_4__0)?
	{ after(grammarAccess.getNodesAccess().getGroup_4()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__Nodes__Group_4__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Nodes__Group_4__0__Impl
	rule__Nodes__Group_4__1
;
finally {
	restoreStackSize(stackSize);
}

rule__Nodes__Group_4__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getNodesAccess().getLeftCurlyBracketKeyword_4_0()); }
	'{'
	{ after(grammarAccess.getNodesAccess().getLeftCurlyBracketKeyword_4_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Nodes__Group_4__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Nodes__Group_4__1__Impl
	rule__Nodes__Group_4__2
;
finally {
	restoreStackSize(stackSize);
}

rule__Nodes__Group_4__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getNodesAccess().getUnorderedGroup_4_1()); }
	(rule__Nodes__UnorderedGroup_4_1)
	{ after(grammarAccess.getNodesAccess().getUnorderedGroup_4_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Nodes__Group_4__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Nodes__Group_4__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__Nodes__Group_4__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getNodesAccess().getRightCurlyBracketKeyword_4_2()); }
	'}'
	{ after(grammarAccess.getNodesAccess().getRightCurlyBracketKeyword_4_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__Size__Group__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Size__Group__0__Impl
	rule__Size__Group__1
;
finally {
	restoreStackSize(stackSize);
}

rule__Size__Group__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getSizeAccess().getSizeAction_0()); }
	()
	{ after(grammarAccess.getSizeAccess().getSizeAction_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Size__Group__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Size__Group__1__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__Size__Group__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getSizeAccess().getGroup_1()); }
	(rule__Size__Group_1__0)
	{ after(grammarAccess.getSizeAccess().getGroup_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__Size__Group_1__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Size__Group_1__0__Impl
	rule__Size__Group_1__1
;
finally {
	restoreStackSize(stackSize);
}

rule__Size__Group_1__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getSizeAccess().getSizeKeyword_1_0()); }
	'size'
	{ after(grammarAccess.getSizeAccess().getSizeKeyword_1_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Size__Group_1__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Size__Group_1__1__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__Size__Group_1__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getSizeAccess().getGroup_1_1()); }
	(rule__Size__Group_1_1__0)?
	{ after(grammarAccess.getSizeAccess().getGroup_1_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__Size__Group_1_1__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Size__Group_1_1__0__Impl
	rule__Size__Group_1_1__1
;
finally {
	restoreStackSize(stackSize);
}

rule__Size__Group_1_1__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getSizeAccess().getLeftCurlyBracketKeyword_1_1_0()); }
	'{'
	{ after(grammarAccess.getSizeAccess().getLeftCurlyBracketKeyword_1_1_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Size__Group_1_1__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Size__Group_1_1__1__Impl
	rule__Size__Group_1_1__2
;
finally {
	restoreStackSize(stackSize);
}

rule__Size__Group_1_1__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getSizeAccess().getUnorderedGroup_1_1_1()); }
	(rule__Size__UnorderedGroup_1_1_1)
	{ after(grammarAccess.getSizeAccess().getUnorderedGroup_1_1_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Size__Group_1_1__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Size__Group_1_1__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__Size__Group_1_1__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getSizeAccess().getRightCurlyBracketKeyword_1_1_2()); }
	'}'
	{ after(grammarAccess.getSizeAccess().getRightCurlyBracketKeyword_1_1_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__Size__Group_1_1_1_0__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Size__Group_1_1_1_0__0__Impl
	rule__Size__Group_1_1_1_0__1
;
finally {
	restoreStackSize(stackSize);
}

rule__Size__Group_1_1_1_0__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getSizeAccess().getHeightKeyword_1_1_1_0_0()); }
	'height'
	{ after(grammarAccess.getSizeAccess().getHeightKeyword_1_1_1_0_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Size__Group_1_1_1_0__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Size__Group_1_1_1_0__1__Impl
	rule__Size__Group_1_1_1_0__2
;
finally {
	restoreStackSize(stackSize);
}

rule__Size__Group_1_1_1_0__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getSizeAccess().getEqualsSignKeyword_1_1_1_0_1()); }
	'='
	{ after(grammarAccess.getSizeAccess().getEqualsSignKeyword_1_1_1_0_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Size__Group_1_1_1_0__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Size__Group_1_1_1_0__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__Size__Group_1_1_1_0__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getSizeAccess().getHeightAssignment_1_1_1_0_2()); }
	(rule__Size__HeightAssignment_1_1_1_0_2)
	{ after(grammarAccess.getSizeAccess().getHeightAssignment_1_1_1_0_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__Size__Group_1_1_1_1__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Size__Group_1_1_1_1__0__Impl
	rule__Size__Group_1_1_1_1__1
;
finally {
	restoreStackSize(stackSize);
}

rule__Size__Group_1_1_1_1__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getSizeAccess().getWidthKeyword_1_1_1_1_0()); }
	'width'
	{ after(grammarAccess.getSizeAccess().getWidthKeyword_1_1_1_1_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Size__Group_1_1_1_1__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Size__Group_1_1_1_1__1__Impl
	rule__Size__Group_1_1_1_1__2
;
finally {
	restoreStackSize(stackSize);
}

rule__Size__Group_1_1_1_1__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getSizeAccess().getEqualsSignKeyword_1_1_1_1_1()); }
	'='
	{ after(grammarAccess.getSizeAccess().getEqualsSignKeyword_1_1_1_1_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Size__Group_1_1_1_1__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Size__Group_1_1_1_1__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__Size__Group_1_1_1_1__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getSizeAccess().getWidthAssignment_1_1_1_1_2()); }
	(rule__Size__WidthAssignment_1_1_1_1_2)
	{ after(grammarAccess.getSizeAccess().getWidthAssignment_1_1_1_1_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__Ports__Group__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Ports__Group__0__Impl
	rule__Ports__Group__1
;
finally {
	restoreStackSize(stackSize);
}

rule__Ports__Group__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getPortsAccess().getPortsAction_0()); }
	()
	{ after(grammarAccess.getPortsAccess().getPortsAction_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Ports__Group__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Ports__Group__1__Impl
	rule__Ports__Group__2
;
finally {
	restoreStackSize(stackSize);
}

rule__Ports__Group__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getPortsAccess().getPortsKeyword_1()); }
	'ports'
	{ after(grammarAccess.getPortsAccess().getPortsKeyword_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Ports__Group__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Ports__Group__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__Ports__Group__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getPortsAccess().getGroup_2()); }
	(rule__Ports__Group_2__0)?
	{ after(grammarAccess.getPortsAccess().getGroup_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__Ports__Group_2__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Ports__Group_2__0__Impl
	rule__Ports__Group_2__1
;
finally {
	restoreStackSize(stackSize);
}

rule__Ports__Group_2__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getPortsAccess().getLeftCurlyBracketKeyword_2_0()); }
	'{'
	{ after(grammarAccess.getPortsAccess().getLeftCurlyBracketKeyword_2_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Ports__Group_2__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Ports__Group_2__1__Impl
	rule__Ports__Group_2__2
;
finally {
	restoreStackSize(stackSize);
}

rule__Ports__Group_2__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getPortsAccess().getUnorderedGroup_2_1()); }
	(rule__Ports__UnorderedGroup_2_1)
	{ after(grammarAccess.getPortsAccess().getUnorderedGroup_2_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Ports__Group_2__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Ports__Group_2__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__Ports__Group_2__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getPortsAccess().getRightCurlyBracketKeyword_2_2()); }
	'}'
	{ after(grammarAccess.getPortsAccess().getRightCurlyBracketKeyword_2_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__Ports__Group_2_1_1__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Ports__Group_2_1_1__0__Impl
	rule__Ports__Group_2_1_1__1
;
finally {
	restoreStackSize(stackSize);
}

rule__Ports__Group_2_1_1__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getPortsAccess().getReUseKeyword_2_1_1_0()); }
	're-use'
	{ after(grammarAccess.getPortsAccess().getReUseKeyword_2_1_1_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Ports__Group_2_1_1__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Ports__Group_2_1_1__1__Impl
	rule__Ports__Group_2_1_1__2
;
finally {
	restoreStackSize(stackSize);
}

rule__Ports__Group_2_1_1__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getPortsAccess().getEqualsSignKeyword_2_1_1_1()); }
	'='
	{ after(grammarAccess.getPortsAccess().getEqualsSignKeyword_2_1_1_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Ports__Group_2_1_1__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Ports__Group_2_1_1__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__Ports__Group_2_1_1__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getPortsAccess().getReUseAssignment_2_1_1_2()); }
	(rule__Ports__ReUseAssignment_2_1_1_2)
	{ after(grammarAccess.getPortsAccess().getReUseAssignment_2_1_1_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__Ports__Group_2_1_3__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Ports__Group_2_1_3__0__Impl
	rule__Ports__Group_2_1_3__1
;
finally {
	restoreStackSize(stackSize);
}

rule__Ports__Group_2_1_3__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getPortsAccess().getConstraintKeyword_2_1_3_0()); }
	'constraint'
	{ after(grammarAccess.getPortsAccess().getConstraintKeyword_2_1_3_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Ports__Group_2_1_3__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Ports__Group_2_1_3__1__Impl
	rule__Ports__Group_2_1_3__2
;
finally {
	restoreStackSize(stackSize);
}

rule__Ports__Group_2_1_3__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getPortsAccess().getEqualsSignKeyword_2_1_3_1()); }
	'='
	{ after(grammarAccess.getPortsAccess().getEqualsSignKeyword_2_1_3_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Ports__Group_2_1_3__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Ports__Group_2_1_3__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__Ports__Group_2_1_3__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getPortsAccess().getConstraintAssignment_2_1_3_2()); }
	(rule__Ports__ConstraintAssignment_2_1_3_2)
	{ after(grammarAccess.getPortsAccess().getConstraintAssignment_2_1_3_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__Flow__Group__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Flow__Group__0__Impl
	rule__Flow__Group__1
;
finally {
	restoreStackSize(stackSize);
}

rule__Flow__Group__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getFlowAccess().getFlowTypeAssignment_0()); }
	(rule__Flow__FlowTypeAssignment_0)
	{ after(grammarAccess.getFlowAccess().getFlowTypeAssignment_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Flow__Group__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Flow__Group__1__Impl
	rule__Flow__Group__2
;
finally {
	restoreStackSize(stackSize);
}

rule__Flow__Group__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getFlowAccess().getSideAssignment_1()); }
	(rule__Flow__SideAssignment_1)
	{ after(grammarAccess.getFlowAccess().getSideAssignment_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Flow__Group__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Flow__Group__2__Impl
	rule__Flow__Group__3
;
finally {
	restoreStackSize(stackSize);
}

rule__Flow__Group__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getFlowAccess().getEqualsSignKeyword_2()); }
	'='
	{ after(grammarAccess.getFlowAccess().getEqualsSignKeyword_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Flow__Group__3
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Flow__Group__3__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__Flow__Group__3__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getFlowAccess().getAmountAssignment_3()); }
	(rule__Flow__AmountAssignment_3)
	{ after(grammarAccess.getFlowAccess().getAmountAssignment_3()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__DoubleQuantity__Group_1__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__DoubleQuantity__Group_1__0__Impl
	rule__DoubleQuantity__Group_1__1
;
finally {
	restoreStackSize(stackSize);
}

rule__DoubleQuantity__Group_1__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getDoubleQuantityAccess().getMinAssignment_1_0()); }
	(rule__DoubleQuantity__MinAssignment_1_0)
	{ after(grammarAccess.getDoubleQuantityAccess().getMinAssignment_1_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__DoubleQuantity__Group_1__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__DoubleQuantity__Group_1__1__Impl
	rule__DoubleQuantity__Group_1__2
;
finally {
	restoreStackSize(stackSize);
}

rule__DoubleQuantity__Group_1__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getDoubleQuantityAccess().getMinMaxAssignment_1_1()); }
	(rule__DoubleQuantity__MinMaxAssignment_1_1)
	{ after(grammarAccess.getDoubleQuantityAccess().getMinMaxAssignment_1_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__DoubleQuantity__Group_1__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__DoubleQuantity__Group_1__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__DoubleQuantity__Group_1__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getDoubleQuantityAccess().getMaxAssignment_1_2()); }
	(rule__DoubleQuantity__MaxAssignment_1_2)
	{ after(grammarAccess.getDoubleQuantityAccess().getMaxAssignment_1_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__DoubleQuantity__Group_2__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__DoubleQuantity__Group_2__0__Impl
	rule__DoubleQuantity__Group_2__1
;
finally {
	restoreStackSize(stackSize);
}

rule__DoubleQuantity__Group_2__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getDoubleQuantityAccess().getMeanAssignment_2_0()); }
	(rule__DoubleQuantity__MeanAssignment_2_0)
	{ after(grammarAccess.getDoubleQuantityAccess().getMeanAssignment_2_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__DoubleQuantity__Group_2__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__DoubleQuantity__Group_2__1__Impl
	rule__DoubleQuantity__Group_2__2
;
finally {
	restoreStackSize(stackSize);
}

rule__DoubleQuantity__Group_2__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getDoubleQuantityAccess().getGaussianAssignment_2_1()); }
	(rule__DoubleQuantity__GaussianAssignment_2_1)
	{ after(grammarAccess.getDoubleQuantityAccess().getGaussianAssignment_2_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__DoubleQuantity__Group_2__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__DoubleQuantity__Group_2__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__DoubleQuantity__Group_2__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getDoubleQuantityAccess().getStddvAssignment_2_2()); }
	(rule__DoubleQuantity__StddvAssignment_2_2)
	{ after(grammarAccess.getDoubleQuantityAccess().getStddvAssignment_2_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__Double__Group__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Double__Group__0__Impl
	rule__Double__Group__1
;
finally {
	restoreStackSize(stackSize);
}

rule__Double__Group__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getDoubleAccess().getINTTerminalRuleCall_0()); }
	RULE_INT
	{ after(grammarAccess.getDoubleAccess().getINTTerminalRuleCall_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Double__Group__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Double__Group__1__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__Double__Group__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getDoubleAccess().getGroup_1()); }
	(rule__Double__Group_1__0)?
	{ after(grammarAccess.getDoubleAccess().getGroup_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__Double__Group_1__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Double__Group_1__0__Impl
	rule__Double__Group_1__1
;
finally {
	restoreStackSize(stackSize);
}

rule__Double__Group_1__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getDoubleAccess().getFullStopKeyword_1_0()); }
	'.'
	{ after(grammarAccess.getDoubleAccess().getFullStopKeyword_1_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Double__Group_1__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Double__Group_1__1__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__Double__Group_1__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getDoubleAccess().getINTTerminalRuleCall_1_1()); }
	RULE_INT
	{ after(grammarAccess.getDoubleAccess().getINTTerminalRuleCall_1_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__Integer__Group__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Integer__Group__0__Impl
	rule__Integer__Group__1
;
finally {
	restoreStackSize(stackSize);
}

rule__Integer__Group__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getIntegerAccess().getINTTerminalRuleCall_0()); }
	RULE_INT
	{ after(grammarAccess.getIntegerAccess().getINTTerminalRuleCall_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Integer__Group__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Integer__Group__1__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__Integer__Group__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getIntegerAccess().getGroup_1()); }
	(rule__Integer__Group_1__0)?
	{ after(grammarAccess.getIntegerAccess().getGroup_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__Integer__Group_1__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Integer__Group_1__0__Impl
	rule__Integer__Group_1__1
;
finally {
	restoreStackSize(stackSize);
}

rule__Integer__Group_1__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getIntegerAccess().getFullStopKeyword_1_0()); }
	'.'
	{ after(grammarAccess.getIntegerAccess().getFullStopKeyword_1_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Integer__Group_1__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Integer__Group_1__1__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__Integer__Group_1__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getIntegerAccess().getINTTerminalRuleCall_1_1()); }
	RULE_INT
	{ after(grammarAccess.getIntegerAccess().getINTTerminalRuleCall_1_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__Configuration__UnorderedGroup_3_1
	@init {
		int stackSize = keepStackSize();
		getUnorderedGroupHelper().enter(grammarAccess.getConfigurationAccess().getUnorderedGroup_3_1());
	}
:
	rule__Configuration__UnorderedGroup_3_1__0
	?
;
finally {
	getUnorderedGroupHelper().leave(grammarAccess.getConfigurationAccess().getUnorderedGroup_3_1());
	restoreStackSize(stackSize);
}

rule__Configuration__UnorderedGroup_3_1__Impl
	@init {
		int stackSize = keepStackSize();
		boolean selected = false;
	}
:
		(
		( 
			{getUnorderedGroupHelper().canSelect(grammarAccess.getConfigurationAccess().getUnorderedGroup_3_1(), 0)}?=>(
				{
					getUnorderedGroupHelper().select(grammarAccess.getConfigurationAccess().getUnorderedGroup_3_1(), 0);
				}
				{
					selected = true;
				}
				(
					{ before(grammarAccess.getConfigurationAccess().getNodesAssignment_3_1_0()); }
					(rule__Configuration__NodesAssignment_3_1_0)
					{ after(grammarAccess.getConfigurationAccess().getNodesAssignment_3_1_0()); }
				)
			)
		)|
		( 
			{getUnorderedGroupHelper().canSelect(grammarAccess.getConfigurationAccess().getUnorderedGroup_3_1(), 1)}?=>(
				{
					getUnorderedGroupHelper().select(grammarAccess.getConfigurationAccess().getUnorderedGroup_3_1(), 1);
				}
				{
					selected = true;
				}
				(
					{ before(grammarAccess.getConfigurationAccess().getEdgesAssignment_3_1_1()); }
					(rule__Configuration__EdgesAssignment_3_1_1)
					{ after(grammarAccess.getConfigurationAccess().getEdgesAssignment_3_1_1()); }
				)
			)
		)|
		( 
			{getUnorderedGroupHelper().canSelect(grammarAccess.getConfigurationAccess().getUnorderedGroup_3_1(), 2)}?=>(
				{
					getUnorderedGroupHelper().select(grammarAccess.getConfigurationAccess().getUnorderedGroup_3_1(), 2);
				}
				{
					selected = true;
				}
				(
					{ before(grammarAccess.getConfigurationAccess().getGroup_3_1_2()); }
					(rule__Configuration__Group_3_1_2__0)
					{ after(grammarAccess.getConfigurationAccess().getGroup_3_1_2()); }
				)
			)
		)|
		( 
			{getUnorderedGroupHelper().canSelect(grammarAccess.getConfigurationAccess().getUnorderedGroup_3_1(), 3)}?=>(
				{
					getUnorderedGroupHelper().select(grammarAccess.getConfigurationAccess().getUnorderedGroup_3_1(), 3);
				}
				{
					selected = true;
				}
				(
					{ before(grammarAccess.getConfigurationAccess().getGroup_3_1_3()); }
					(rule__Configuration__Group_3_1_3__0)
					{ after(grammarAccess.getConfigurationAccess().getGroup_3_1_3()); }
				)
			)
		)|
		( 
			{getUnorderedGroupHelper().canSelect(grammarAccess.getConfigurationAccess().getUnorderedGroup_3_1(), 4)}?=>(
				{
					getUnorderedGroupHelper().select(grammarAccess.getConfigurationAccess().getUnorderedGroup_3_1(), 4);
				}
				{
					selected = true;
				}
				(
					{ before(grammarAccess.getConfigurationAccess().getGroup_3_1_4()); }
					(rule__Configuration__Group_3_1_4__0)
					{ after(grammarAccess.getConfigurationAccess().getGroup_3_1_4()); }
				)
			)
		)|
		( 
			{getUnorderedGroupHelper().canSelect(grammarAccess.getConfigurationAccess().getUnorderedGroup_3_1(), 5)}?=>(
				{
					getUnorderedGroupHelper().select(grammarAccess.getConfigurationAccess().getUnorderedGroup_3_1(), 5);
				}
				{
					selected = true;
				}
				(
					{ before(grammarAccess.getConfigurationAccess().getHierarchyAssignment_3_1_5()); }
					(rule__Configuration__HierarchyAssignment_3_1_5)
					{ after(grammarAccess.getConfigurationAccess().getHierarchyAssignment_3_1_5()); }
				)
			)
		)|
		( 
			{getUnorderedGroupHelper().canSelect(grammarAccess.getConfigurationAccess().getUnorderedGroup_3_1(), 6)}?=>(
				{
					getUnorderedGroupHelper().select(grammarAccess.getConfigurationAccess().getUnorderedGroup_3_1(), 6);
				}
				{
					selected = true;
				}
				(
					{ before(grammarAccess.getConfigurationAccess().getGroup_3_1_6()); }
					(rule__Configuration__Group_3_1_6__0)
					{ after(grammarAccess.getConfigurationAccess().getGroup_3_1_6()); }
				)
			)
		)|
		( 
			{getUnorderedGroupHelper().canSelect(grammarAccess.getConfigurationAccess().getUnorderedGroup_3_1(), 7)}?=>(
				{
					getUnorderedGroupHelper().select(grammarAccess.getConfigurationAccess().getUnorderedGroup_3_1(), 7);
				}
				{
					selected = true;
				}
				(
					{ before(grammarAccess.getConfigurationAccess().getGroup_3_1_7()); }
					(rule__Configuration__Group_3_1_7__0)
					{ after(grammarAccess.getConfigurationAccess().getGroup_3_1_7()); }
				)
			)
		)|
		( 
			{getUnorderedGroupHelper().canSelect(grammarAccess.getConfigurationAccess().getUnorderedGroup_3_1(), 8)}?=>(
				{
					getUnorderedGroupHelper().select(grammarAccess.getConfigurationAccess().getUnorderedGroup_3_1(), 8);
				}
				{
					selected = true;
				}
				(
					{ before(grammarAccess.getConfigurationAccess().getGroup_3_1_8()); }
					(rule__Configuration__Group_3_1_8__0)
					{ after(grammarAccess.getConfigurationAccess().getGroup_3_1_8()); }
				)
			)
		)
		)
;
finally {
	if (selected)
		getUnorderedGroupHelper().returnFromSelection(grammarAccess.getConfigurationAccess().getUnorderedGroup_3_1());
	restoreStackSize(stackSize);
}

rule__Configuration__UnorderedGroup_3_1__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Configuration__UnorderedGroup_3_1__Impl
	rule__Configuration__UnorderedGroup_3_1__1?
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__UnorderedGroup_3_1__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Configuration__UnorderedGroup_3_1__Impl
	rule__Configuration__UnorderedGroup_3_1__2?
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__UnorderedGroup_3_1__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Configuration__UnorderedGroup_3_1__Impl
	rule__Configuration__UnorderedGroup_3_1__3?
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__UnorderedGroup_3_1__3
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Configuration__UnorderedGroup_3_1__Impl
	rule__Configuration__UnorderedGroup_3_1__4?
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__UnorderedGroup_3_1__4
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Configuration__UnorderedGroup_3_1__Impl
	rule__Configuration__UnorderedGroup_3_1__5?
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__UnorderedGroup_3_1__5
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Configuration__UnorderedGroup_3_1__Impl
	rule__Configuration__UnorderedGroup_3_1__6?
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__UnorderedGroup_3_1__6
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Configuration__UnorderedGroup_3_1__Impl
	rule__Configuration__UnorderedGroup_3_1__7?
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__UnorderedGroup_3_1__7
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Configuration__UnorderedGroup_3_1__Impl
	rule__Configuration__UnorderedGroup_3_1__8?
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__UnorderedGroup_3_1__8
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Configuration__UnorderedGroup_3_1__Impl
;
finally {
	restoreStackSize(stackSize);
}


rule__Hierarchy__UnorderedGroup_2_1
	@init {
		int stackSize = keepStackSize();
		getUnorderedGroupHelper().enter(grammarAccess.getHierarchyAccess().getUnorderedGroup_2_1());
	}
:
	rule__Hierarchy__UnorderedGroup_2_1__0
	?
;
finally {
	getUnorderedGroupHelper().leave(grammarAccess.getHierarchyAccess().getUnorderedGroup_2_1());
	restoreStackSize(stackSize);
}

rule__Hierarchy__UnorderedGroup_2_1__Impl
	@init {
		int stackSize = keepStackSize();
		boolean selected = false;
	}
:
		(
		( 
			{getUnorderedGroupHelper().canSelect(grammarAccess.getHierarchyAccess().getUnorderedGroup_2_1(), 0)}?=>(
				{
					getUnorderedGroupHelper().select(grammarAccess.getHierarchyAccess().getUnorderedGroup_2_1(), 0);
				}
				{
					selected = true;
				}
				(
					{ before(grammarAccess.getHierarchyAccess().getGroup_2_1_0()); }
					(rule__Hierarchy__Group_2_1_0__0)
					{ after(grammarAccess.getHierarchyAccess().getGroup_2_1_0()); }
				)
			)
		)|
		( 
			{getUnorderedGroupHelper().canSelect(grammarAccess.getHierarchyAccess().getUnorderedGroup_2_1(), 1)}?=>(
				{
					getUnorderedGroupHelper().select(grammarAccess.getHierarchyAccess().getUnorderedGroup_2_1(), 1);
				}
				{
					selected = true;
				}
				(
					{ before(grammarAccess.getHierarchyAccess().getGroup_2_1_1()); }
					(rule__Hierarchy__Group_2_1_1__0)
					{ after(grammarAccess.getHierarchyAccess().getGroup_2_1_1()); }
				)
			)
		)|
		( 
			{getUnorderedGroupHelper().canSelect(grammarAccess.getHierarchyAccess().getUnorderedGroup_2_1(), 2)}?=>(
				{
					getUnorderedGroupHelper().select(grammarAccess.getHierarchyAccess().getUnorderedGroup_2_1(), 2);
				}
				{
					selected = true;
				}
				(
					{ before(grammarAccess.getHierarchyAccess().getGroup_2_1_2()); }
					(rule__Hierarchy__Group_2_1_2__0)
					{ after(grammarAccess.getHierarchyAccess().getGroup_2_1_2()); }
				)
			)
		)|
		( 
			{getUnorderedGroupHelper().canSelect(grammarAccess.getHierarchyAccess().getUnorderedGroup_2_1(), 3)}?=>(
				{
					getUnorderedGroupHelper().select(grammarAccess.getHierarchyAccess().getUnorderedGroup_2_1(), 3);
				}
				{
					selected = true;
				}
				(
					{ before(grammarAccess.getHierarchyAccess().getGroup_2_1_3()); }
					(rule__Hierarchy__Group_2_1_3__0)
					{ after(grammarAccess.getHierarchyAccess().getGroup_2_1_3()); }
				)
			)
		)
		)
;
finally {
	if (selected)
		getUnorderedGroupHelper().returnFromSelection(grammarAccess.getHierarchyAccess().getUnorderedGroup_2_1());
	restoreStackSize(stackSize);
}

rule__Hierarchy__UnorderedGroup_2_1__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Hierarchy__UnorderedGroup_2_1__Impl
	rule__Hierarchy__UnorderedGroup_2_1__1?
;
finally {
	restoreStackSize(stackSize);
}

rule__Hierarchy__UnorderedGroup_2_1__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Hierarchy__UnorderedGroup_2_1__Impl
	rule__Hierarchy__UnorderedGroup_2_1__2?
;
finally {
	restoreStackSize(stackSize);
}

rule__Hierarchy__UnorderedGroup_2_1__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Hierarchy__UnorderedGroup_2_1__Impl
	rule__Hierarchy__UnorderedGroup_2_1__3?
;
finally {
	restoreStackSize(stackSize);
}

rule__Hierarchy__UnorderedGroup_2_1__3
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Hierarchy__UnorderedGroup_2_1__Impl
;
finally {
	restoreStackSize(stackSize);
}


rule__Edges__UnorderedGroup_1_1
	@init {
		int stackSize = keepStackSize();
		getUnorderedGroupHelper().enter(grammarAccess.getEdgesAccess().getUnorderedGroup_1_1());
	}
:
	rule__Edges__UnorderedGroup_1_1__0
	?
;
finally {
	getUnorderedGroupHelper().leave(grammarAccess.getEdgesAccess().getUnorderedGroup_1_1());
	restoreStackSize(stackSize);
}

rule__Edges__UnorderedGroup_1_1__Impl
	@init {
		int stackSize = keepStackSize();
		boolean selected = false;
	}
:
		(
		( 
			{getUnorderedGroupHelper().canSelect(grammarAccess.getEdgesAccess().getUnorderedGroup_1_1(), 0)}?=>(
				{
					getUnorderedGroupHelper().select(grammarAccess.getEdgesAccess().getUnorderedGroup_1_1(), 0);
				}
				{
					selected = true;
				}
				(
					{ before(grammarAccess.getEdgesAccess().getLabelsAssignment_1_1_0()); }
					(rule__Edges__LabelsAssignment_1_1_0)
					{ after(grammarAccess.getEdgesAccess().getLabelsAssignment_1_1_0()); }
				)
			)
		)|
		( 
			{getUnorderedGroupHelper().canSelect(grammarAccess.getEdgesAccess().getUnorderedGroup_1_1(), 1)}?=>(
				{
					getUnorderedGroupHelper().select(grammarAccess.getEdgesAccess().getUnorderedGroup_1_1(), 1);
				}
				{
					selected = true;
				}
				(
					{ before(grammarAccess.getEdgesAccess().getSelfLoopsAssignment_1_1_1()); }
					(rule__Edges__SelfLoopsAssignment_1_1_1)
					{ after(grammarAccess.getEdgesAccess().getSelfLoopsAssignment_1_1_1()); }
				)
			)
		)
		)
;
finally {
	if (selected)
		getUnorderedGroupHelper().returnFromSelection(grammarAccess.getEdgesAccess().getUnorderedGroup_1_1());
	restoreStackSize(stackSize);
}

rule__Edges__UnorderedGroup_1_1__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Edges__UnorderedGroup_1_1__Impl
	rule__Edges__UnorderedGroup_1_1__1?
;
finally {
	restoreStackSize(stackSize);
}

rule__Edges__UnorderedGroup_1_1__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Edges__UnorderedGroup_1_1__Impl
;
finally {
	restoreStackSize(stackSize);
}


rule__Nodes__UnorderedGroup_4_1
	@init {
		int stackSize = keepStackSize();
		getUnorderedGroupHelper().enter(grammarAccess.getNodesAccess().getUnorderedGroup_4_1());
	}
:
	rule__Nodes__UnorderedGroup_4_1__0
	?
;
finally {
	getUnorderedGroupHelper().leave(grammarAccess.getNodesAccess().getUnorderedGroup_4_1());
	restoreStackSize(stackSize);
}

rule__Nodes__UnorderedGroup_4_1__Impl
	@init {
		int stackSize = keepStackSize();
		boolean selected = false;
	}
:
		(
		( 
			{getUnorderedGroupHelper().canSelect(grammarAccess.getNodesAccess().getUnorderedGroup_4_1(), 0)}?=>(
				{
					getUnorderedGroupHelper().select(grammarAccess.getNodesAccess().getUnorderedGroup_4_1(), 0);
				}
				{
					selected = true;
				}
				(
					{ before(grammarAccess.getNodesAccess().getPortsAssignment_4_1_0()); }
					(rule__Nodes__PortsAssignment_4_1_0)
					{ after(grammarAccess.getNodesAccess().getPortsAssignment_4_1_0()); }
				)
			)
		)|
		( 
			{getUnorderedGroupHelper().canSelect(grammarAccess.getNodesAccess().getUnorderedGroup_4_1(), 1)}?=>(
				{
					getUnorderedGroupHelper().select(grammarAccess.getNodesAccess().getUnorderedGroup_4_1(), 1);
				}
				{
					selected = true;
				}
				(
					{ before(grammarAccess.getNodesAccess().getLabelsAssignment_4_1_1()); }
					(rule__Nodes__LabelsAssignment_4_1_1)
					{ after(grammarAccess.getNodesAccess().getLabelsAssignment_4_1_1()); }
				)
			)
		)|
		( 
			{getUnorderedGroupHelper().canSelect(grammarAccess.getNodesAccess().getUnorderedGroup_4_1(), 2)}?=>(
				{
					getUnorderedGroupHelper().select(grammarAccess.getNodesAccess().getUnorderedGroup_4_1(), 2);
				}
				{
					selected = true;
				}
				(
					{ before(grammarAccess.getNodesAccess().getSizeAssignment_4_1_2()); }
					(rule__Nodes__SizeAssignment_4_1_2)
					{ after(grammarAccess.getNodesAccess().getSizeAssignment_4_1_2()); }
				)
			)
		)|
		( 
			{getUnorderedGroupHelper().canSelect(grammarAccess.getNodesAccess().getUnorderedGroup_4_1(), 3)}?=>(
				{
					getUnorderedGroupHelper().select(grammarAccess.getNodesAccess().getUnorderedGroup_4_1(), 3);
				}
				{
					selected = true;
				}
				(
					{ before(grammarAccess.getNodesAccess().getRemoveIsolatedAssignment_4_1_3()); }
					(rule__Nodes__RemoveIsolatedAssignment_4_1_3)
					{ after(grammarAccess.getNodesAccess().getRemoveIsolatedAssignment_4_1_3()); }
				)
			)
		)
		)
;
finally {
	if (selected)
		getUnorderedGroupHelper().returnFromSelection(grammarAccess.getNodesAccess().getUnorderedGroup_4_1());
	restoreStackSize(stackSize);
}

rule__Nodes__UnorderedGroup_4_1__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Nodes__UnorderedGroup_4_1__Impl
	rule__Nodes__UnorderedGroup_4_1__1?
;
finally {
	restoreStackSize(stackSize);
}

rule__Nodes__UnorderedGroup_4_1__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Nodes__UnorderedGroup_4_1__Impl
	rule__Nodes__UnorderedGroup_4_1__2?
;
finally {
	restoreStackSize(stackSize);
}

rule__Nodes__UnorderedGroup_4_1__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Nodes__UnorderedGroup_4_1__Impl
	rule__Nodes__UnorderedGroup_4_1__3?
;
finally {
	restoreStackSize(stackSize);
}

rule__Nodes__UnorderedGroup_4_1__3
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Nodes__UnorderedGroup_4_1__Impl
;
finally {
	restoreStackSize(stackSize);
}


rule__Size__UnorderedGroup_1_1_1
	@init {
		int stackSize = keepStackSize();
		getUnorderedGroupHelper().enter(grammarAccess.getSizeAccess().getUnorderedGroup_1_1_1());
	}
:
	rule__Size__UnorderedGroup_1_1_1__0
	?
;
finally {
	getUnorderedGroupHelper().leave(grammarAccess.getSizeAccess().getUnorderedGroup_1_1_1());
	restoreStackSize(stackSize);
}

rule__Size__UnorderedGroup_1_1_1__Impl
	@init {
		int stackSize = keepStackSize();
		boolean selected = false;
	}
:
		(
		( 
			{getUnorderedGroupHelper().canSelect(grammarAccess.getSizeAccess().getUnorderedGroup_1_1_1(), 0)}?=>(
				{
					getUnorderedGroupHelper().select(grammarAccess.getSizeAccess().getUnorderedGroup_1_1_1(), 0);
				}
				{
					selected = true;
				}
				(
					{ before(grammarAccess.getSizeAccess().getGroup_1_1_1_0()); }
					(rule__Size__Group_1_1_1_0__0)
					{ after(grammarAccess.getSizeAccess().getGroup_1_1_1_0()); }
				)
			)
		)|
		( 
			{getUnorderedGroupHelper().canSelect(grammarAccess.getSizeAccess().getUnorderedGroup_1_1_1(), 1)}?=>(
				{
					getUnorderedGroupHelper().select(grammarAccess.getSizeAccess().getUnorderedGroup_1_1_1(), 1);
				}
				{
					selected = true;
				}
				(
					{ before(grammarAccess.getSizeAccess().getGroup_1_1_1_1()); }
					(rule__Size__Group_1_1_1_1__0)
					{ after(grammarAccess.getSizeAccess().getGroup_1_1_1_1()); }
				)
			)
		)
		)
;
finally {
	if (selected)
		getUnorderedGroupHelper().returnFromSelection(grammarAccess.getSizeAccess().getUnorderedGroup_1_1_1());
	restoreStackSize(stackSize);
}

rule__Size__UnorderedGroup_1_1_1__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Size__UnorderedGroup_1_1_1__Impl
	rule__Size__UnorderedGroup_1_1_1__1?
;
finally {
	restoreStackSize(stackSize);
}

rule__Size__UnorderedGroup_1_1_1__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Size__UnorderedGroup_1_1_1__Impl
;
finally {
	restoreStackSize(stackSize);
}


rule__Ports__UnorderedGroup_2_1
	@init {
		int stackSize = keepStackSize();
		getUnorderedGroupHelper().enter(grammarAccess.getPortsAccess().getUnorderedGroup_2_1());
	}
:
	rule__Ports__UnorderedGroup_2_1__0
	?
;
finally {
	getUnorderedGroupHelper().leave(grammarAccess.getPortsAccess().getUnorderedGroup_2_1());
	restoreStackSize(stackSize);
}

rule__Ports__UnorderedGroup_2_1__Impl
	@init {
		int stackSize = keepStackSize();
		boolean selected = false;
	}
:
		(
		( 
			{getUnorderedGroupHelper().canSelect(grammarAccess.getPortsAccess().getUnorderedGroup_2_1(), 0)}?=>(
				{
					getUnorderedGroupHelper().select(grammarAccess.getPortsAccess().getUnorderedGroup_2_1(), 0);
				}
				{
					selected = true;
				}
				(
					{ before(grammarAccess.getPortsAccess().getLabelsAssignment_2_1_0()); }
					(rule__Ports__LabelsAssignment_2_1_0)
					{ after(grammarAccess.getPortsAccess().getLabelsAssignment_2_1_0()); }
				)
			)
		)|
		( 
			{getUnorderedGroupHelper().canSelect(grammarAccess.getPortsAccess().getUnorderedGroup_2_1(), 1)}?=>(
				{
					getUnorderedGroupHelper().select(grammarAccess.getPortsAccess().getUnorderedGroup_2_1(), 1);
				}
				{
					selected = true;
				}
				(
					{ before(grammarAccess.getPortsAccess().getGroup_2_1_1()); }
					(rule__Ports__Group_2_1_1__0)
					{ after(grammarAccess.getPortsAccess().getGroup_2_1_1()); }
				)
			)
		)|
		( 
			{getUnorderedGroupHelper().canSelect(grammarAccess.getPortsAccess().getUnorderedGroup_2_1(), 2)}?=>(
				{
					getUnorderedGroupHelper().select(grammarAccess.getPortsAccess().getUnorderedGroup_2_1(), 2);
				}
				{
					selected = true;
				}
				(
					{ before(grammarAccess.getPortsAccess().getSizeAssignment_2_1_2()); }
					(rule__Ports__SizeAssignment_2_1_2)
					{ after(grammarAccess.getPortsAccess().getSizeAssignment_2_1_2()); }
				)
			)
		)|
		( 
			{getUnorderedGroupHelper().canSelect(grammarAccess.getPortsAccess().getUnorderedGroup_2_1(), 3)}?=>(
				{
					getUnorderedGroupHelper().select(grammarAccess.getPortsAccess().getUnorderedGroup_2_1(), 3);
				}
				{
					selected = true;
				}
				(
					{ before(grammarAccess.getPortsAccess().getGroup_2_1_3()); }
					(rule__Ports__Group_2_1_3__0)
					{ after(grammarAccess.getPortsAccess().getGroup_2_1_3()); }
				)
			)
		)|
		( 
			{getUnorderedGroupHelper().canSelect(grammarAccess.getPortsAccess().getUnorderedGroup_2_1(), 4)}?=>(
				{
					getUnorderedGroupHelper().select(grammarAccess.getPortsAccess().getUnorderedGroup_2_1(), 4);
				}
				{
					selected = true;
				}
				(
					(
						{ before(grammarAccess.getPortsAccess().getFlowAssignment_2_1_4()); }
						(rule__Ports__FlowAssignment_2_1_4)
						{ after(grammarAccess.getPortsAccess().getFlowAssignment_2_1_4()); }
					)
					(
						{ before(grammarAccess.getPortsAccess().getFlowAssignment_2_1_4()); }
						((rule__Ports__FlowAssignment_2_1_4)=>rule__Ports__FlowAssignment_2_1_4)*
						{ after(grammarAccess.getPortsAccess().getFlowAssignment_2_1_4()); }
					)
				)
			)
		)
		)
;
finally {
	if (selected)
		getUnorderedGroupHelper().returnFromSelection(grammarAccess.getPortsAccess().getUnorderedGroup_2_1());
	restoreStackSize(stackSize);
}

rule__Ports__UnorderedGroup_2_1__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Ports__UnorderedGroup_2_1__Impl
	rule__Ports__UnorderedGroup_2_1__1?
;
finally {
	restoreStackSize(stackSize);
}

rule__Ports__UnorderedGroup_2_1__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Ports__UnorderedGroup_2_1__Impl
	rule__Ports__UnorderedGroup_2_1__2?
;
finally {
	restoreStackSize(stackSize);
}

rule__Ports__UnorderedGroup_2_1__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Ports__UnorderedGroup_2_1__Impl
	rule__Ports__UnorderedGroup_2_1__3?
;
finally {
	restoreStackSize(stackSize);
}

rule__Ports__UnorderedGroup_2_1__3
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Ports__UnorderedGroup_2_1__Impl
	rule__Ports__UnorderedGroup_2_1__4?
;
finally {
	restoreStackSize(stackSize);
}

rule__Ports__UnorderedGroup_2_1__4
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Ports__UnorderedGroup_2_1__Impl
;
finally {
	restoreStackSize(stackSize);
}


rule__RandGraph__ConfigsAssignment
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getRandGraphAccess().getConfigsConfigurationParserRuleCall_0()); }
		ruleConfiguration
		{ after(grammarAccess.getRandGraphAccess().getConfigsConfigurationParserRuleCall_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__SamplesAssignment_1
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getConfigurationAccess().getSamplesINTTerminalRuleCall_1_0()); }
		RULE_INT
		{ after(grammarAccess.getConfigurationAccess().getSamplesINTTerminalRuleCall_1_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__FormAssignment_2
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getConfigurationAccess().getFormFormEnumRuleCall_2_0()); }
		ruleForm
		{ after(grammarAccess.getConfigurationAccess().getFormFormEnumRuleCall_2_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__NodesAssignment_3_1_0
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getConfigurationAccess().getNodesNodesParserRuleCall_3_1_0_0()); }
		ruleNodes
		{ after(grammarAccess.getConfigurationAccess().getNodesNodesParserRuleCall_3_1_0_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__EdgesAssignment_3_1_1
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getConfigurationAccess().getEdgesEdgesParserRuleCall_3_1_1_0()); }
		ruleEdges
		{ after(grammarAccess.getConfigurationAccess().getEdgesEdgesParserRuleCall_3_1_1_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__MWAssignment_3_1_2_0
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getConfigurationAccess().getMWMaxWidthKeyword_3_1_2_0_0()); }
		(
			{ before(grammarAccess.getConfigurationAccess().getMWMaxWidthKeyword_3_1_2_0_0()); }
			'maxWidth'
			{ after(grammarAccess.getConfigurationAccess().getMWMaxWidthKeyword_3_1_2_0_0()); }
		)
		{ after(grammarAccess.getConfigurationAccess().getMWMaxWidthKeyword_3_1_2_0_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__MaxWidthAssignment_3_1_2_2
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getConfigurationAccess().getMaxWidthIntegerParserRuleCall_3_1_2_2_0()); }
		ruleInteger
		{ after(grammarAccess.getConfigurationAccess().getMaxWidthIntegerParserRuleCall_3_1_2_2_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__MDAssignment_3_1_3_0
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getConfigurationAccess().getMDMaxDegreeKeyword_3_1_3_0_0()); }
		(
			{ before(grammarAccess.getConfigurationAccess().getMDMaxDegreeKeyword_3_1_3_0_0()); }
			'maxDegree'
			{ after(grammarAccess.getConfigurationAccess().getMDMaxDegreeKeyword_3_1_3_0_0()); }
		)
		{ after(grammarAccess.getConfigurationAccess().getMDMaxDegreeKeyword_3_1_3_0_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__MaxDegreeAssignment_3_1_3_2
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getConfigurationAccess().getMaxDegreeIntegerParserRuleCall_3_1_3_2_0()); }
		ruleInteger
		{ after(grammarAccess.getConfigurationAccess().getMaxDegreeIntegerParserRuleCall_3_1_3_2_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__PFAssignment_3_1_4_0
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getConfigurationAccess().getPFPartitionFractionKeyword_3_1_4_0_0()); }
		(
			{ before(grammarAccess.getConfigurationAccess().getPFPartitionFractionKeyword_3_1_4_0_0()); }
			'partitionFraction'
			{ after(grammarAccess.getConfigurationAccess().getPFPartitionFractionKeyword_3_1_4_0_0()); }
		)
		{ after(grammarAccess.getConfigurationAccess().getPFPartitionFractionKeyword_3_1_4_0_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__FractionAssignment_3_1_4_2
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getConfigurationAccess().getFractionDoubleQuantityParserRuleCall_3_1_4_2_0()); }
		ruleDoubleQuantity
		{ after(grammarAccess.getConfigurationAccess().getFractionDoubleQuantityParserRuleCall_3_1_4_2_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__HierarchyAssignment_3_1_5
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getConfigurationAccess().getHierarchyHierarchyParserRuleCall_3_1_5_0()); }
		ruleHierarchy
		{ after(grammarAccess.getConfigurationAccess().getHierarchyHierarchyParserRuleCall_3_1_5_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__SeedAssignment_3_1_6_2
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getConfigurationAccess().getSeedIntegerParserRuleCall_3_1_6_2_0()); }
		ruleInteger
		{ after(grammarAccess.getConfigurationAccess().getSeedIntegerParserRuleCall_3_1_6_2_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__FormatAssignment_3_1_7_2
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getConfigurationAccess().getFormatFormatsEnumRuleCall_3_1_7_2_0()); }
		ruleFormats
		{ after(grammarAccess.getConfigurationAccess().getFormatFormatsEnumRuleCall_3_1_7_2_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__Configuration__FilenameAssignment_3_1_8_2
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getConfigurationAccess().getFilenameSTRINGTerminalRuleCall_3_1_8_2_0()); }
		RULE_STRING
		{ after(grammarAccess.getConfigurationAccess().getFilenameSTRINGTerminalRuleCall_3_1_8_2_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__Hierarchy__LevelsAssignment_2_1_0_2
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getHierarchyAccess().getLevelsDoubleQuantityParserRuleCall_2_1_0_2_0()); }
		ruleDoubleQuantity
		{ after(grammarAccess.getHierarchyAccess().getLevelsDoubleQuantityParserRuleCall_2_1_0_2_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__Hierarchy__EdgesAssignment_2_1_1_2
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getHierarchyAccess().getEdgesDoubleQuantityParserRuleCall_2_1_1_2_0()); }
		ruleDoubleQuantity
		{ after(grammarAccess.getHierarchyAccess().getEdgesDoubleQuantityParserRuleCall_2_1_1_2_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__Hierarchy__NumHierarchNodesAssignment_2_1_2_2
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getHierarchyAccess().getNumHierarchNodesDoubleQuantityParserRuleCall_2_1_2_2_0()); }
		ruleDoubleQuantity
		{ after(grammarAccess.getHierarchyAccess().getNumHierarchNodesDoubleQuantityParserRuleCall_2_1_2_2_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__Hierarchy__CrossHierarchRelAssignment_2_1_3_2
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getHierarchyAccess().getCrossHierarchRelDoubleQuantityParserRuleCall_2_1_3_2_0()); }
		ruleDoubleQuantity
		{ after(grammarAccess.getHierarchyAccess().getCrossHierarchRelDoubleQuantityParserRuleCall_2_1_3_2_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__Edges__DensityAssignment_0_1_0
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getEdgesAccess().getDensityDensityKeyword_0_1_0_0()); }
		(
			{ before(grammarAccess.getEdgesAccess().getDensityDensityKeyword_0_1_0_0()); }
			'density'
			{ after(grammarAccess.getEdgesAccess().getDensityDensityKeyword_0_1_0_0()); }
		)
		{ after(grammarAccess.getEdgesAccess().getDensityDensityKeyword_0_1_0_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__Edges__TotalAssignment_0_1_1
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getEdgesAccess().getTotalTotalKeyword_0_1_1_0()); }
		(
			{ before(grammarAccess.getEdgesAccess().getTotalTotalKeyword_0_1_1_0()); }
			'total'
			{ after(grammarAccess.getEdgesAccess().getTotalTotalKeyword_0_1_1_0()); }
		)
		{ after(grammarAccess.getEdgesAccess().getTotalTotalKeyword_0_1_1_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__Edges__RelativeAssignment_0_1_2
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getEdgesAccess().getRelativeRelativeKeyword_0_1_2_0()); }
		(
			{ before(grammarAccess.getEdgesAccess().getRelativeRelativeKeyword_0_1_2_0()); }
			'relative'
			{ after(grammarAccess.getEdgesAccess().getRelativeRelativeKeyword_0_1_2_0()); }
		)
		{ after(grammarAccess.getEdgesAccess().getRelativeRelativeKeyword_0_1_2_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__Edges__OutboundAssignment_0_1_3
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getEdgesAccess().getOutboundOutgoingKeyword_0_1_3_0()); }
		(
			{ before(grammarAccess.getEdgesAccess().getOutboundOutgoingKeyword_0_1_3_0()); }
			'outgoing'
			{ after(grammarAccess.getEdgesAccess().getOutboundOutgoingKeyword_0_1_3_0()); }
		)
		{ after(grammarAccess.getEdgesAccess().getOutboundOutgoingKeyword_0_1_3_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__Edges__NEdgesAssignment_0_3
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getEdgesAccess().getNEdgesDoubleQuantityParserRuleCall_0_3_0()); }
		ruleDoubleQuantity
		{ after(grammarAccess.getEdgesAccess().getNEdgesDoubleQuantityParserRuleCall_0_3_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__Edges__LabelsAssignment_1_1_0
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getEdgesAccess().getLabelsLabelsKeyword_1_1_0_0()); }
		(
			{ before(grammarAccess.getEdgesAccess().getLabelsLabelsKeyword_1_1_0_0()); }
			'labels'
			{ after(grammarAccess.getEdgesAccess().getLabelsLabelsKeyword_1_1_0_0()); }
		)
		{ after(grammarAccess.getEdgesAccess().getLabelsLabelsKeyword_1_1_0_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__Edges__SelfLoopsAssignment_1_1_1
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getEdgesAccess().getSelfLoopsSelfLoopsKeyword_1_1_1_0()); }
		(
			{ before(grammarAccess.getEdgesAccess().getSelfLoopsSelfLoopsKeyword_1_1_1_0()); }
			'self loops'
			{ after(grammarAccess.getEdgesAccess().getSelfLoopsSelfLoopsKeyword_1_1_1_0()); }
		)
		{ after(grammarAccess.getEdgesAccess().getSelfLoopsSelfLoopsKeyword_1_1_1_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__Nodes__NNodesAssignment_3
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getNodesAccess().getNNodesDoubleQuantityParserRuleCall_3_0()); }
		ruleDoubleQuantity
		{ after(grammarAccess.getNodesAccess().getNNodesDoubleQuantityParserRuleCall_3_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__Nodes__PortsAssignment_4_1_0
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getNodesAccess().getPortsPortsParserRuleCall_4_1_0_0()); }
		rulePorts
		{ after(grammarAccess.getNodesAccess().getPortsPortsParserRuleCall_4_1_0_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__Nodes__LabelsAssignment_4_1_1
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getNodesAccess().getLabelsLabelsParserRuleCall_4_1_1_0()); }
		ruleLabels
		{ after(grammarAccess.getNodesAccess().getLabelsLabelsParserRuleCall_4_1_1_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__Nodes__SizeAssignment_4_1_2
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getNodesAccess().getSizeSizeParserRuleCall_4_1_2_0()); }
		ruleSize
		{ after(grammarAccess.getNodesAccess().getSizeSizeParserRuleCall_4_1_2_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__Nodes__RemoveIsolatedAssignment_4_1_3
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getNodesAccess().getRemoveIsolatedRemoveIsolatedKeyword_4_1_3_0()); }
		(
			{ before(grammarAccess.getNodesAccess().getRemoveIsolatedRemoveIsolatedKeyword_4_1_3_0()); }
			'remove isolated'
			{ after(grammarAccess.getNodesAccess().getRemoveIsolatedRemoveIsolatedKeyword_4_1_3_0()); }
		)
		{ after(grammarAccess.getNodesAccess().getRemoveIsolatedRemoveIsolatedKeyword_4_1_3_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__Size__HeightAssignment_1_1_1_0_2
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getSizeAccess().getHeightDoubleQuantityParserRuleCall_1_1_1_0_2_0()); }
		ruleDoubleQuantity
		{ after(grammarAccess.getSizeAccess().getHeightDoubleQuantityParserRuleCall_1_1_1_0_2_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__Size__WidthAssignment_1_1_1_1_2
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getSizeAccess().getWidthDoubleQuantityParserRuleCall_1_1_1_1_2_0()); }
		ruleDoubleQuantity
		{ after(grammarAccess.getSizeAccess().getWidthDoubleQuantityParserRuleCall_1_1_1_1_2_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__Ports__LabelsAssignment_2_1_0
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getPortsAccess().getLabelsLabelsParserRuleCall_2_1_0_0()); }
		ruleLabels
		{ after(grammarAccess.getPortsAccess().getLabelsLabelsParserRuleCall_2_1_0_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__Ports__ReUseAssignment_2_1_1_2
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getPortsAccess().getReUseDoubleQuantityParserRuleCall_2_1_1_2_0()); }
		ruleDoubleQuantity
		{ after(grammarAccess.getPortsAccess().getReUseDoubleQuantityParserRuleCall_2_1_1_2_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__Ports__SizeAssignment_2_1_2
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getPortsAccess().getSizeSizeParserRuleCall_2_1_2_0()); }
		ruleSize
		{ after(grammarAccess.getPortsAccess().getSizeSizeParserRuleCall_2_1_2_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__Ports__ConstraintAssignment_2_1_3_2
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getPortsAccess().getConstraintConstraintTypeEnumRuleCall_2_1_3_2_0()); }
		ruleConstraintType
		{ after(grammarAccess.getPortsAccess().getConstraintConstraintTypeEnumRuleCall_2_1_3_2_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__Ports__FlowAssignment_2_1_4
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getPortsAccess().getFlowFlowParserRuleCall_2_1_4_0()); }
		ruleFlow
		{ after(grammarAccess.getPortsAccess().getFlowFlowParserRuleCall_2_1_4_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__Flow__FlowTypeAssignment_0
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getFlowAccess().getFlowTypeFlowTypeEnumRuleCall_0_0()); }
		ruleFlowType
		{ after(grammarAccess.getFlowAccess().getFlowTypeFlowTypeEnumRuleCall_0_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__Flow__SideAssignment_1
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getFlowAccess().getSideSideEnumRuleCall_1_0()); }
		ruleSide
		{ after(grammarAccess.getFlowAccess().getSideSideEnumRuleCall_1_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__Flow__AmountAssignment_3
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getFlowAccess().getAmountDoubleQuantityParserRuleCall_3_0()); }
		ruleDoubleQuantity
		{ after(grammarAccess.getFlowAccess().getAmountDoubleQuantityParserRuleCall_3_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__DoubleQuantity__QuantAssignment_0
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getDoubleQuantityAccess().getQuantDoubleParserRuleCall_0_0()); }
		ruleDouble
		{ after(grammarAccess.getDoubleQuantityAccess().getQuantDoubleParserRuleCall_0_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__DoubleQuantity__MinAssignment_1_0
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getDoubleQuantityAccess().getMinDoubleParserRuleCall_1_0_0()); }
		ruleDouble
		{ after(grammarAccess.getDoubleQuantityAccess().getMinDoubleParserRuleCall_1_0_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__DoubleQuantity__MinMaxAssignment_1_1
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getDoubleQuantityAccess().getMinMaxToKeyword_1_1_0()); }
		(
			{ before(grammarAccess.getDoubleQuantityAccess().getMinMaxToKeyword_1_1_0()); }
			'to'
			{ after(grammarAccess.getDoubleQuantityAccess().getMinMaxToKeyword_1_1_0()); }
		)
		{ after(grammarAccess.getDoubleQuantityAccess().getMinMaxToKeyword_1_1_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__DoubleQuantity__MaxAssignment_1_2
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getDoubleQuantityAccess().getMaxDoubleParserRuleCall_1_2_0()); }
		ruleDouble
		{ after(grammarAccess.getDoubleQuantityAccess().getMaxDoubleParserRuleCall_1_2_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__DoubleQuantity__MeanAssignment_2_0
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getDoubleQuantityAccess().getMeanDoubleParserRuleCall_2_0_0()); }
		ruleDouble
		{ after(grammarAccess.getDoubleQuantityAccess().getMeanDoubleParserRuleCall_2_0_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__DoubleQuantity__GaussianAssignment_2_1
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getDoubleQuantityAccess().getGaussianPmParserRuleCall_2_1_0()); }
		rulePm
		{ after(grammarAccess.getDoubleQuantityAccess().getGaussianPmParserRuleCall_2_1_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__DoubleQuantity__StddvAssignment_2_2
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getDoubleQuantityAccess().getStddvDoubleParserRuleCall_2_2_0()); }
		ruleDouble
		{ after(grammarAccess.getDoubleQuantityAccess().getStddvDoubleParserRuleCall_2_2_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

RULE_ID : '^'? ('a'..'z'|'A'..'Z'|'_') ('a'..'z'|'A'..'Z'|'_'|'0'..'9')*;

RULE_INT : ('0'..'9')+;

RULE_STRING : ('"' ('\\' .|~(('\\'|'"')))* '"'|'\'' ('\\' .|~(('\\'|'\'')))* '\'');

RULE_ML_COMMENT : '/*' ( options {greedy=false;} : . )*'*/';

RULE_SL_COMMENT : '//' ~(('\n'|'\r'))* ('\r'? '\n')?;

RULE_WS : (' '|'\t'|'\r'|'\n')+;

RULE_ANY_OTHER : .;
