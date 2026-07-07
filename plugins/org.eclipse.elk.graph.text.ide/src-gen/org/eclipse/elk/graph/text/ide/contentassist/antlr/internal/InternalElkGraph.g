/*******************************************************************************
 * Copyright (c) 2016 TypeFox GmbH (http://www.typefox.io) and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
grammar InternalElkGraph;

options {
	superClass=AbstractInternalContentAssistParser;
}

@lexer::header {
package org.eclipse.elk.graph.text.ide.contentassist.antlr.internal;

// Hack: Use our own Lexer superclass by means of import. 
// Currently there is no other way to specify the superclass for the lexer.
import org.eclipse.xtext.ide.editor.contentassist.antlr.internal.Lexer;
}

@parser::header {
package org.eclipse.elk.graph.text.ide.contentassist.antlr.internal;

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
import org.eclipse.elk.graph.text.services.ElkGraphGrammarAccess;

}
@parser::members {
	private ElkGraphGrammarAccess grammarAccess;

	public void setGrammarAccess(ElkGraphGrammarAccess grammarAccess) {
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

// Entry rule entryRuleRootNode
entryRuleRootNode
:
{ before(grammarAccess.getRootNodeRule()); }
	 ruleRootNode
{ after(grammarAccess.getRootNodeRule()); } 
	 EOF 
;

// Rule RootNode
ruleRootNode 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getRootNodeAccess().getGroup()); }
		(rule__RootNode__Group__0)
		{ after(grammarAccess.getRootNodeAccess().getGroup()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

// Entry rule entryRuleElkNode
entryRuleElkNode
:
{ before(grammarAccess.getElkNodeRule()); }
	 ruleElkNode
{ after(grammarAccess.getElkNodeRule()); } 
	 EOF 
;

// Rule ElkNode
ruleElkNode 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getElkNodeAccess().getGroup()); }
		(rule__ElkNode__Group__0)
		{ after(grammarAccess.getElkNodeAccess().getGroup()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

// Entry rule entryRuleElkLabel
entryRuleElkLabel
:
{ before(grammarAccess.getElkLabelRule()); }
	 ruleElkLabel
{ after(grammarAccess.getElkLabelRule()); } 
	 EOF 
;

// Rule ElkLabel
ruleElkLabel 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getElkLabelAccess().getGroup()); }
		(rule__ElkLabel__Group__0)
		{ after(grammarAccess.getElkLabelAccess().getGroup()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

// Entry rule entryRuleElkPort
entryRuleElkPort
:
{ before(grammarAccess.getElkPortRule()); }
	 ruleElkPort
{ after(grammarAccess.getElkPortRule()); } 
	 EOF 
;

// Rule ElkPort
ruleElkPort 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getElkPortAccess().getGroup()); }
		(rule__ElkPort__Group__0)
		{ after(grammarAccess.getElkPortAccess().getGroup()); }
	)
;
finally {
	restoreStackSize(stackSize);
}


// Rule ShapeLayout
ruleShapeLayout 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getShapeLayoutAccess().getGroup()); }
		(rule__ShapeLayout__Group__0)
		{ after(grammarAccess.getShapeLayoutAccess().getGroup()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

// Entry rule entryRuleElkEdge
entryRuleElkEdge
:
{ before(grammarAccess.getElkEdgeRule()); }
	 ruleElkEdge
{ after(grammarAccess.getElkEdgeRule()); } 
	 EOF 
;

// Rule ElkEdge
ruleElkEdge 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getElkEdgeAccess().getGroup()); }
		(rule__ElkEdge__Group__0)
		{ after(grammarAccess.getElkEdgeAccess().getGroup()); }
	)
;
finally {
	restoreStackSize(stackSize);
}


// Rule EdgeLayout
ruleEdgeLayout 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getEdgeLayoutAccess().getGroup()); }
		(rule__EdgeLayout__Group__0)
		{ after(grammarAccess.getEdgeLayoutAccess().getGroup()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

// Entry rule entryRuleElkSingleEdgeSection
entryRuleElkSingleEdgeSection
:
{ before(grammarAccess.getElkSingleEdgeSectionRule()); }
	 ruleElkSingleEdgeSection
{ after(grammarAccess.getElkSingleEdgeSectionRule()); } 
	 EOF 
;

// Rule ElkSingleEdgeSection
ruleElkSingleEdgeSection 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getElkSingleEdgeSectionAccess().getGroup()); }
		(rule__ElkSingleEdgeSection__Group__0)
		{ after(grammarAccess.getElkSingleEdgeSectionAccess().getGroup()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

// Entry rule entryRuleElkEdgeSection
entryRuleElkEdgeSection
:
{ before(grammarAccess.getElkEdgeSectionRule()); }
	 ruleElkEdgeSection
{ after(grammarAccess.getElkEdgeSectionRule()); } 
	 EOF 
;

// Rule ElkEdgeSection
ruleElkEdgeSection 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getElkEdgeSectionAccess().getGroup()); }
		(rule__ElkEdgeSection__Group__0)
		{ after(grammarAccess.getElkEdgeSectionAccess().getGroup()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

// Entry rule entryRuleElkBendPoint
entryRuleElkBendPoint
:
{ before(grammarAccess.getElkBendPointRule()); }
	 ruleElkBendPoint
{ after(grammarAccess.getElkBendPointRule()); } 
	 EOF 
;

// Rule ElkBendPoint
ruleElkBendPoint 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getElkBendPointAccess().getGroup()); }
		(rule__ElkBendPoint__Group__0)
		{ after(grammarAccess.getElkBendPointAccess().getGroup()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

// Entry rule entryRuleQualifiedId
entryRuleQualifiedId
:
{ before(grammarAccess.getQualifiedIdRule()); }
	 ruleQualifiedId
{ after(grammarAccess.getQualifiedIdRule()); } 
	 EOF 
;

// Rule QualifiedId
ruleQualifiedId 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getQualifiedIdAccess().getGroup()); }
		(rule__QualifiedId__Group__0)
		{ after(grammarAccess.getQualifiedIdAccess().getGroup()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

// Entry rule entryRuleNumber
entryRuleNumber
:
{ before(grammarAccess.getNumberRule()); }
	 ruleNumber
{ after(grammarAccess.getNumberRule()); } 
	 EOF 
;

// Rule Number
ruleNumber 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getNumberAccess().getAlternatives()); }
		(rule__Number__Alternatives)
		{ after(grammarAccess.getNumberAccess().getAlternatives()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

// Entry rule entryRuleProperty
entryRuleProperty
:
{ before(grammarAccess.getPropertyRule()); }
	 ruleProperty
{ after(grammarAccess.getPropertyRule()); } 
	 EOF 
;

// Rule Property
ruleProperty 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getPropertyAccess().getGroup()); }
		(rule__Property__Group__0)
		{ after(grammarAccess.getPropertyAccess().getGroup()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

// Entry rule entryRulePropertyKey
entryRulePropertyKey
@init { 
	HiddenTokens myHiddenTokenState = ((XtextTokenStream)input).setHiddenTokens();
}
:
{ before(grammarAccess.getPropertyKeyRule()); }
	 rulePropertyKey
{ after(grammarAccess.getPropertyKeyRule()); } 
	 EOF 
;
finally {
	myHiddenTokenState.restore();
}

// Rule PropertyKey
rulePropertyKey 
	@init {
		HiddenTokens myHiddenTokenState = ((XtextTokenStream)input).setHiddenTokens();
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getPropertyKeyAccess().getGroup()); }
		(rule__PropertyKey__Group__0)
		{ after(grammarAccess.getPropertyKeyAccess().getGroup()); }
	)
;
finally {
	restoreStackSize(stackSize);
	myHiddenTokenState.restore();
}

// Entry rule entryRuleStringValue
entryRuleStringValue
:
{ before(grammarAccess.getStringValueRule()); }
	 ruleStringValue
{ after(grammarAccess.getStringValueRule()); } 
	 EOF 
;

// Rule StringValue
ruleStringValue 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getStringValueAccess().getSTRINGTerminalRuleCall()); }
		RULE_STRING
		{ after(grammarAccess.getStringValueAccess().getSTRINGTerminalRuleCall()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

// Entry rule entryRuleQualifiedIdValue
entryRuleQualifiedIdValue
:
{ before(grammarAccess.getQualifiedIdValueRule()); }
	 ruleQualifiedIdValue
{ after(grammarAccess.getQualifiedIdValueRule()); } 
	 EOF 
;

// Rule QualifiedIdValue
ruleQualifiedIdValue 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getQualifiedIdValueAccess().getQualifiedIdParserRuleCall()); }
		ruleQualifiedId
		{ after(grammarAccess.getQualifiedIdValueAccess().getQualifiedIdParserRuleCall()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

// Entry rule entryRuleNumberValue
entryRuleNumberValue
:
{ before(grammarAccess.getNumberValueRule()); }
	 ruleNumberValue
{ after(grammarAccess.getNumberValueRule()); } 
	 EOF 
;

// Rule NumberValue
ruleNumberValue 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getNumberValueAccess().getAlternatives()); }
		(rule__NumberValue__Alternatives)
		{ after(grammarAccess.getNumberValueAccess().getAlternatives()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

// Entry rule entryRuleBooleanValue
entryRuleBooleanValue
:
{ before(grammarAccess.getBooleanValueRule()); }
	 ruleBooleanValue
{ after(grammarAccess.getBooleanValueRule()); } 
	 EOF 
;

// Rule BooleanValue
ruleBooleanValue 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getBooleanValueAccess().getAlternatives()); }
		(rule__BooleanValue__Alternatives)
		{ after(grammarAccess.getBooleanValueAccess().getAlternatives()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__RootNode__Alternatives_4
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getRootNodeAccess().getLabelsAssignment_4_0()); }
		(rule__RootNode__LabelsAssignment_4_0)
		{ after(grammarAccess.getRootNodeAccess().getLabelsAssignment_4_0()); }
	)
	|
	(
		{ before(grammarAccess.getRootNodeAccess().getPortsAssignment_4_1()); }
		(rule__RootNode__PortsAssignment_4_1)
		{ after(grammarAccess.getRootNodeAccess().getPortsAssignment_4_1()); }
	)
	|
	(
		{ before(grammarAccess.getRootNodeAccess().getChildrenAssignment_4_2()); }
		(rule__RootNode__ChildrenAssignment_4_2)
		{ after(grammarAccess.getRootNodeAccess().getChildrenAssignment_4_2()); }
	)
	|
	(
		{ before(grammarAccess.getRootNodeAccess().getContainedEdgesAssignment_4_3()); }
		(rule__RootNode__ContainedEdgesAssignment_4_3)
		{ after(grammarAccess.getRootNodeAccess().getContainedEdgesAssignment_4_3()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNode__Alternatives_2_3
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkNodeAccess().getLabelsAssignment_2_3_0()); }
		(rule__ElkNode__LabelsAssignment_2_3_0)
		{ after(grammarAccess.getElkNodeAccess().getLabelsAssignment_2_3_0()); }
	)
	|
	(
		{ before(grammarAccess.getElkNodeAccess().getPortsAssignment_2_3_1()); }
		(rule__ElkNode__PortsAssignment_2_3_1)
		{ after(grammarAccess.getElkNodeAccess().getPortsAssignment_2_3_1()); }
	)
	|
	(
		{ before(grammarAccess.getElkNodeAccess().getChildrenAssignment_2_3_2()); }
		(rule__ElkNode__ChildrenAssignment_2_3_2)
		{ after(grammarAccess.getElkNodeAccess().getChildrenAssignment_2_3_2()); }
	)
	|
	(
		{ before(grammarAccess.getElkNodeAccess().getContainedEdgesAssignment_2_3_3()); }
		(rule__ElkNode__ContainedEdgesAssignment_2_3_3)
		{ after(grammarAccess.getElkNodeAccess().getContainedEdgesAssignment_2_3_3()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__EdgeLayout__Alternatives_2
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getEdgeLayoutAccess().getSectionsAssignment_2_0()); }
		(rule__EdgeLayout__SectionsAssignment_2_0)
		{ after(grammarAccess.getEdgeLayoutAccess().getSectionsAssignment_2_0()); }
	)
	|
	(
		(
			{ before(grammarAccess.getEdgeLayoutAccess().getSectionsAssignment_2_1()); }
			(rule__EdgeLayout__SectionsAssignment_2_1)
			{ after(grammarAccess.getEdgeLayoutAccess().getSectionsAssignment_2_1()); }
		)
		(
			{ before(grammarAccess.getEdgeLayoutAccess().getSectionsAssignment_2_1()); }
			(rule__EdgeLayout__SectionsAssignment_2_1)*
			{ after(grammarAccess.getEdgeLayoutAccess().getSectionsAssignment_2_1()); }
		)
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__Number__Alternatives
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getNumberAccess().getSIGNED_INTTerminalRuleCall_0()); }
		RULE_SIGNED_INT
		{ after(grammarAccess.getNumberAccess().getSIGNED_INTTerminalRuleCall_0()); }
	)
	|
	(
		{ before(grammarAccess.getNumberAccess().getFLOATTerminalRuleCall_1()); }
		RULE_FLOAT
		{ after(grammarAccess.getNumberAccess().getFLOATTerminalRuleCall_1()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__Property__Alternatives_2
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getPropertyAccess().getValueAssignment_2_0()); }
		(rule__Property__ValueAssignment_2_0)
		{ after(grammarAccess.getPropertyAccess().getValueAssignment_2_0()); }
	)
	|
	(
		{ before(grammarAccess.getPropertyAccess().getValueAssignment_2_1()); }
		(rule__Property__ValueAssignment_2_1)
		{ after(grammarAccess.getPropertyAccess().getValueAssignment_2_1()); }
	)
	|
	(
		{ before(grammarAccess.getPropertyAccess().getValueAssignment_2_2()); }
		(rule__Property__ValueAssignment_2_2)
		{ after(grammarAccess.getPropertyAccess().getValueAssignment_2_2()); }
	)
	|
	(
		{ before(grammarAccess.getPropertyAccess().getValueAssignment_2_3()); }
		(rule__Property__ValueAssignment_2_3)
		{ after(grammarAccess.getPropertyAccess().getValueAssignment_2_3()); }
	)
	|
	(
		{ before(grammarAccess.getPropertyAccess().getNullKeyword_2_4()); }
		'null'
		{ after(grammarAccess.getPropertyAccess().getNullKeyword_2_4()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__NumberValue__Alternatives
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getNumberValueAccess().getSIGNED_INTTerminalRuleCall_0()); }
		RULE_SIGNED_INT
		{ after(grammarAccess.getNumberValueAccess().getSIGNED_INTTerminalRuleCall_0()); }
	)
	|
	(
		{ before(grammarAccess.getNumberValueAccess().getFLOATTerminalRuleCall_1()); }
		RULE_FLOAT
		{ after(grammarAccess.getNumberValueAccess().getFLOATTerminalRuleCall_1()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__BooleanValue__Alternatives
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getBooleanValueAccess().getTrueKeyword_0()); }
		'true'
		{ after(grammarAccess.getBooleanValueAccess().getTrueKeyword_0()); }
	)
	|
	(
		{ before(grammarAccess.getBooleanValueAccess().getFalseKeyword_1()); }
		'false'
		{ after(grammarAccess.getBooleanValueAccess().getFalseKeyword_1()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__RootNode__Group__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__RootNode__Group__0__Impl
	rule__RootNode__Group__1
;
finally {
	restoreStackSize(stackSize);
}

rule__RootNode__Group__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getRootNodeAccess().getElkNodeAction_0()); }
	()
	{ after(grammarAccess.getRootNodeAccess().getElkNodeAction_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__RootNode__Group__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__RootNode__Group__1__Impl
	rule__RootNode__Group__2
;
finally {
	restoreStackSize(stackSize);
}

rule__RootNode__Group__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getRootNodeAccess().getGroup_1()); }
	(rule__RootNode__Group_1__0)?
	{ after(grammarAccess.getRootNodeAccess().getGroup_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__RootNode__Group__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__RootNode__Group__2__Impl
	rule__RootNode__Group__3
;
finally {
	restoreStackSize(stackSize);
}

rule__RootNode__Group__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getRootNodeAccess().getShapeLayoutParserRuleCall_2()); }
	(ruleShapeLayout)?
	{ after(grammarAccess.getRootNodeAccess().getShapeLayoutParserRuleCall_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__RootNode__Group__3
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__RootNode__Group__3__Impl
	rule__RootNode__Group__4
;
finally {
	restoreStackSize(stackSize);
}

rule__RootNode__Group__3__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getRootNodeAccess().getPropertiesAssignment_3()); }
	(rule__RootNode__PropertiesAssignment_3)*
	{ after(grammarAccess.getRootNodeAccess().getPropertiesAssignment_3()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__RootNode__Group__4
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__RootNode__Group__4__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__RootNode__Group__4__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getRootNodeAccess().getAlternatives_4()); }
	(rule__RootNode__Alternatives_4)*
	{ after(grammarAccess.getRootNodeAccess().getAlternatives_4()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__RootNode__Group_1__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__RootNode__Group_1__0__Impl
	rule__RootNode__Group_1__1
;
finally {
	restoreStackSize(stackSize);
}

rule__RootNode__Group_1__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getRootNodeAccess().getGraphKeyword_1_0()); }
	'graph'
	{ after(grammarAccess.getRootNodeAccess().getGraphKeyword_1_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__RootNode__Group_1__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__RootNode__Group_1__1__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__RootNode__Group_1__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getRootNodeAccess().getIdentifierAssignment_1_1()); }
	(rule__RootNode__IdentifierAssignment_1_1)
	{ after(grammarAccess.getRootNodeAccess().getIdentifierAssignment_1_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkNode__Group__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkNode__Group__0__Impl
	rule__ElkNode__Group__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNode__Group__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkNodeAccess().getNodeKeyword_0()); }
	'node'
	{ after(grammarAccess.getElkNodeAccess().getNodeKeyword_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNode__Group__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkNode__Group__1__Impl
	rule__ElkNode__Group__2
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNode__Group__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkNodeAccess().getIdentifierAssignment_1()); }
	(rule__ElkNode__IdentifierAssignment_1)
	{ after(grammarAccess.getElkNodeAccess().getIdentifierAssignment_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNode__Group__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkNode__Group__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNode__Group__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkNodeAccess().getGroup_2()); }
	(rule__ElkNode__Group_2__0)?
	{ after(grammarAccess.getElkNodeAccess().getGroup_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkNode__Group_2__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkNode__Group_2__0__Impl
	rule__ElkNode__Group_2__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNode__Group_2__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkNodeAccess().getLeftCurlyBracketKeyword_2_0()); }
	'{'
	{ after(grammarAccess.getElkNodeAccess().getLeftCurlyBracketKeyword_2_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNode__Group_2__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkNode__Group_2__1__Impl
	rule__ElkNode__Group_2__2
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNode__Group_2__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkNodeAccess().getShapeLayoutParserRuleCall_2_1()); }
	(ruleShapeLayout)?
	{ after(grammarAccess.getElkNodeAccess().getShapeLayoutParserRuleCall_2_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNode__Group_2__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkNode__Group_2__2__Impl
	rule__ElkNode__Group_2__3
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNode__Group_2__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkNodeAccess().getPropertiesAssignment_2_2()); }
	(rule__ElkNode__PropertiesAssignment_2_2)*
	{ after(grammarAccess.getElkNodeAccess().getPropertiesAssignment_2_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNode__Group_2__3
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkNode__Group_2__3__Impl
	rule__ElkNode__Group_2__4
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNode__Group_2__3__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkNodeAccess().getAlternatives_2_3()); }
	(rule__ElkNode__Alternatives_2_3)*
	{ after(grammarAccess.getElkNodeAccess().getAlternatives_2_3()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNode__Group_2__4
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkNode__Group_2__4__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNode__Group_2__4__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkNodeAccess().getRightCurlyBracketKeyword_2_4()); }
	'}'
	{ after(grammarAccess.getElkNodeAccess().getRightCurlyBracketKeyword_2_4()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkLabel__Group__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkLabel__Group__0__Impl
	rule__ElkLabel__Group__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkLabel__Group__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkLabelAccess().getLabelKeyword_0()); }
	'label'
	{ after(grammarAccess.getElkLabelAccess().getLabelKeyword_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkLabel__Group__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkLabel__Group__1__Impl
	rule__ElkLabel__Group__2
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkLabel__Group__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkLabelAccess().getGroup_1()); }
	(rule__ElkLabel__Group_1__0)?
	{ after(grammarAccess.getElkLabelAccess().getGroup_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkLabel__Group__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkLabel__Group__2__Impl
	rule__ElkLabel__Group__3
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkLabel__Group__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkLabelAccess().getTextAssignment_2()); }
	(rule__ElkLabel__TextAssignment_2)
	{ after(grammarAccess.getElkLabelAccess().getTextAssignment_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkLabel__Group__3
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkLabel__Group__3__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkLabel__Group__3__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkLabelAccess().getGroup_3()); }
	(rule__ElkLabel__Group_3__0)?
	{ after(grammarAccess.getElkLabelAccess().getGroup_3()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkLabel__Group_1__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkLabel__Group_1__0__Impl
	rule__ElkLabel__Group_1__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkLabel__Group_1__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkLabelAccess().getIdentifierAssignment_1_0()); }
	(rule__ElkLabel__IdentifierAssignment_1_0)
	{ after(grammarAccess.getElkLabelAccess().getIdentifierAssignment_1_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkLabel__Group_1__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkLabel__Group_1__1__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkLabel__Group_1__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkLabelAccess().getColonKeyword_1_1()); }
	':'
	{ after(grammarAccess.getElkLabelAccess().getColonKeyword_1_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkLabel__Group_3__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkLabel__Group_3__0__Impl
	rule__ElkLabel__Group_3__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkLabel__Group_3__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkLabelAccess().getLeftCurlyBracketKeyword_3_0()); }
	'{'
	{ after(grammarAccess.getElkLabelAccess().getLeftCurlyBracketKeyword_3_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkLabel__Group_3__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkLabel__Group_3__1__Impl
	rule__ElkLabel__Group_3__2
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkLabel__Group_3__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkLabelAccess().getShapeLayoutParserRuleCall_3_1()); }
	(ruleShapeLayout)?
	{ after(grammarAccess.getElkLabelAccess().getShapeLayoutParserRuleCall_3_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkLabel__Group_3__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkLabel__Group_3__2__Impl
	rule__ElkLabel__Group_3__3
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkLabel__Group_3__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkLabelAccess().getPropertiesAssignment_3_2()); }
	(rule__ElkLabel__PropertiesAssignment_3_2)*
	{ after(grammarAccess.getElkLabelAccess().getPropertiesAssignment_3_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkLabel__Group_3__3
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkLabel__Group_3__3__Impl
	rule__ElkLabel__Group_3__4
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkLabel__Group_3__3__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkLabelAccess().getLabelsAssignment_3_3()); }
	(rule__ElkLabel__LabelsAssignment_3_3)*
	{ after(grammarAccess.getElkLabelAccess().getLabelsAssignment_3_3()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkLabel__Group_3__4
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkLabel__Group_3__4__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkLabel__Group_3__4__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkLabelAccess().getRightCurlyBracketKeyword_3_4()); }
	'}'
	{ after(grammarAccess.getElkLabelAccess().getRightCurlyBracketKeyword_3_4()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkPort__Group__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkPort__Group__0__Impl
	rule__ElkPort__Group__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkPort__Group__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkPortAccess().getPortKeyword_0()); }
	'port'
	{ after(grammarAccess.getElkPortAccess().getPortKeyword_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkPort__Group__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkPort__Group__1__Impl
	rule__ElkPort__Group__2
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkPort__Group__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkPortAccess().getIdentifierAssignment_1()); }
	(rule__ElkPort__IdentifierAssignment_1)
	{ after(grammarAccess.getElkPortAccess().getIdentifierAssignment_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkPort__Group__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkPort__Group__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkPort__Group__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkPortAccess().getGroup_2()); }
	(rule__ElkPort__Group_2__0)?
	{ after(grammarAccess.getElkPortAccess().getGroup_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkPort__Group_2__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkPort__Group_2__0__Impl
	rule__ElkPort__Group_2__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkPort__Group_2__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkPortAccess().getLeftCurlyBracketKeyword_2_0()); }
	'{'
	{ after(grammarAccess.getElkPortAccess().getLeftCurlyBracketKeyword_2_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkPort__Group_2__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkPort__Group_2__1__Impl
	rule__ElkPort__Group_2__2
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkPort__Group_2__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkPortAccess().getShapeLayoutParserRuleCall_2_1()); }
	(ruleShapeLayout)?
	{ after(grammarAccess.getElkPortAccess().getShapeLayoutParserRuleCall_2_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkPort__Group_2__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkPort__Group_2__2__Impl
	rule__ElkPort__Group_2__3
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkPort__Group_2__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkPortAccess().getPropertiesAssignment_2_2()); }
	(rule__ElkPort__PropertiesAssignment_2_2)*
	{ after(grammarAccess.getElkPortAccess().getPropertiesAssignment_2_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkPort__Group_2__3
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkPort__Group_2__3__Impl
	rule__ElkPort__Group_2__4
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkPort__Group_2__3__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkPortAccess().getLabelsAssignment_2_3()); }
	(rule__ElkPort__LabelsAssignment_2_3)*
	{ after(grammarAccess.getElkPortAccess().getLabelsAssignment_2_3()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkPort__Group_2__4
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkPort__Group_2__4__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkPort__Group_2__4__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkPortAccess().getRightCurlyBracketKeyword_2_4()); }
	'}'
	{ after(grammarAccess.getElkPortAccess().getRightCurlyBracketKeyword_2_4()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ShapeLayout__Group__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ShapeLayout__Group__0__Impl
	rule__ShapeLayout__Group__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeLayout__Group__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getShapeLayoutAccess().getLayoutKeyword_0()); }
	'layout'
	{ after(grammarAccess.getShapeLayoutAccess().getLayoutKeyword_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeLayout__Group__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ShapeLayout__Group__1__Impl
	rule__ShapeLayout__Group__2
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeLayout__Group__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getShapeLayoutAccess().getLeftSquareBracketKeyword_1()); }
	'['
	{ after(grammarAccess.getShapeLayoutAccess().getLeftSquareBracketKeyword_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeLayout__Group__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ShapeLayout__Group__2__Impl
	rule__ShapeLayout__Group__3
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeLayout__Group__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getShapeLayoutAccess().getUnorderedGroup_2()); }
	(rule__ShapeLayout__UnorderedGroup_2)
	{ after(grammarAccess.getShapeLayoutAccess().getUnorderedGroup_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeLayout__Group__3
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ShapeLayout__Group__3__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeLayout__Group__3__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getShapeLayoutAccess().getRightSquareBracketKeyword_3()); }
	']'
	{ after(grammarAccess.getShapeLayoutAccess().getRightSquareBracketKeyword_3()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ShapeLayout__Group_2_0__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ShapeLayout__Group_2_0__0__Impl
	rule__ShapeLayout__Group_2_0__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeLayout__Group_2_0__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getShapeLayoutAccess().getPositionKeyword_2_0_0()); }
	'position'
	{ after(grammarAccess.getShapeLayoutAccess().getPositionKeyword_2_0_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeLayout__Group_2_0__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ShapeLayout__Group_2_0__1__Impl
	rule__ShapeLayout__Group_2_0__2
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeLayout__Group_2_0__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getShapeLayoutAccess().getColonKeyword_2_0_1()); }
	':'
	{ after(grammarAccess.getShapeLayoutAccess().getColonKeyword_2_0_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeLayout__Group_2_0__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ShapeLayout__Group_2_0__2__Impl
	rule__ShapeLayout__Group_2_0__3
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeLayout__Group_2_0__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getShapeLayoutAccess().getXAssignment_2_0_2()); }
	(rule__ShapeLayout__XAssignment_2_0_2)
	{ after(grammarAccess.getShapeLayoutAccess().getXAssignment_2_0_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeLayout__Group_2_0__3
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ShapeLayout__Group_2_0__3__Impl
	rule__ShapeLayout__Group_2_0__4
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeLayout__Group_2_0__3__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getShapeLayoutAccess().getCommaKeyword_2_0_3()); }
	','
	{ after(grammarAccess.getShapeLayoutAccess().getCommaKeyword_2_0_3()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeLayout__Group_2_0__4
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ShapeLayout__Group_2_0__4__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeLayout__Group_2_0__4__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getShapeLayoutAccess().getYAssignment_2_0_4()); }
	(rule__ShapeLayout__YAssignment_2_0_4)
	{ after(grammarAccess.getShapeLayoutAccess().getYAssignment_2_0_4()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ShapeLayout__Group_2_1__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ShapeLayout__Group_2_1__0__Impl
	rule__ShapeLayout__Group_2_1__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeLayout__Group_2_1__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getShapeLayoutAccess().getSizeKeyword_2_1_0()); }
	'size'
	{ after(grammarAccess.getShapeLayoutAccess().getSizeKeyword_2_1_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeLayout__Group_2_1__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ShapeLayout__Group_2_1__1__Impl
	rule__ShapeLayout__Group_2_1__2
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeLayout__Group_2_1__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getShapeLayoutAccess().getColonKeyword_2_1_1()); }
	':'
	{ after(grammarAccess.getShapeLayoutAccess().getColonKeyword_2_1_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeLayout__Group_2_1__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ShapeLayout__Group_2_1__2__Impl
	rule__ShapeLayout__Group_2_1__3
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeLayout__Group_2_1__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getShapeLayoutAccess().getWidthAssignment_2_1_2()); }
	(rule__ShapeLayout__WidthAssignment_2_1_2)
	{ after(grammarAccess.getShapeLayoutAccess().getWidthAssignment_2_1_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeLayout__Group_2_1__3
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ShapeLayout__Group_2_1__3__Impl
	rule__ShapeLayout__Group_2_1__4
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeLayout__Group_2_1__3__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getShapeLayoutAccess().getCommaKeyword_2_1_3()); }
	','
	{ after(grammarAccess.getShapeLayoutAccess().getCommaKeyword_2_1_3()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeLayout__Group_2_1__4
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ShapeLayout__Group_2_1__4__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeLayout__Group_2_1__4__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getShapeLayoutAccess().getHeightAssignment_2_1_4()); }
	(rule__ShapeLayout__HeightAssignment_2_1_4)
	{ after(grammarAccess.getShapeLayoutAccess().getHeightAssignment_2_1_4()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkEdge__Group__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdge__Group__0__Impl
	rule__ElkEdge__Group__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdge__Group__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeAccess().getEdgeKeyword_0()); }
	'edge'
	{ after(grammarAccess.getElkEdgeAccess().getEdgeKeyword_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdge__Group__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdge__Group__1__Impl
	rule__ElkEdge__Group__2
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdge__Group__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeAccess().getGroup_1()); }
	(rule__ElkEdge__Group_1__0)?
	{ after(grammarAccess.getElkEdgeAccess().getGroup_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdge__Group__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdge__Group__2__Impl
	rule__ElkEdge__Group__3
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdge__Group__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeAccess().getSourcesAssignment_2()); }
	(rule__ElkEdge__SourcesAssignment_2)
	{ after(grammarAccess.getElkEdgeAccess().getSourcesAssignment_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdge__Group__3
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdge__Group__3__Impl
	rule__ElkEdge__Group__4
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdge__Group__3__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeAccess().getGroup_3()); }
	(rule__ElkEdge__Group_3__0)*
	{ after(grammarAccess.getElkEdgeAccess().getGroup_3()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdge__Group__4
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdge__Group__4__Impl
	rule__ElkEdge__Group__5
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdge__Group__4__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeAccess().getHyphenMinusGreaterThanSignKeyword_4()); }
	'->'
	{ after(grammarAccess.getElkEdgeAccess().getHyphenMinusGreaterThanSignKeyword_4()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdge__Group__5
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdge__Group__5__Impl
	rule__ElkEdge__Group__6
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdge__Group__5__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeAccess().getTargetsAssignment_5()); }
	(rule__ElkEdge__TargetsAssignment_5)
	{ after(grammarAccess.getElkEdgeAccess().getTargetsAssignment_5()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdge__Group__6
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdge__Group__6__Impl
	rule__ElkEdge__Group__7
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdge__Group__6__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeAccess().getGroup_6()); }
	(rule__ElkEdge__Group_6__0)*
	{ after(grammarAccess.getElkEdgeAccess().getGroup_6()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdge__Group__7
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdge__Group__7__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdge__Group__7__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeAccess().getGroup_7()); }
	(rule__ElkEdge__Group_7__0)?
	{ after(grammarAccess.getElkEdgeAccess().getGroup_7()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkEdge__Group_1__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdge__Group_1__0__Impl
	rule__ElkEdge__Group_1__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdge__Group_1__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeAccess().getIdentifierAssignment_1_0()); }
	(rule__ElkEdge__IdentifierAssignment_1_0)
	{ after(grammarAccess.getElkEdgeAccess().getIdentifierAssignment_1_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdge__Group_1__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdge__Group_1__1__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdge__Group_1__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeAccess().getColonKeyword_1_1()); }
	':'
	{ after(grammarAccess.getElkEdgeAccess().getColonKeyword_1_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkEdge__Group_3__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdge__Group_3__0__Impl
	rule__ElkEdge__Group_3__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdge__Group_3__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeAccess().getCommaKeyword_3_0()); }
	','
	{ after(grammarAccess.getElkEdgeAccess().getCommaKeyword_3_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdge__Group_3__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdge__Group_3__1__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdge__Group_3__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeAccess().getSourcesAssignment_3_1()); }
	(rule__ElkEdge__SourcesAssignment_3_1)
	{ after(grammarAccess.getElkEdgeAccess().getSourcesAssignment_3_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkEdge__Group_6__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdge__Group_6__0__Impl
	rule__ElkEdge__Group_6__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdge__Group_6__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeAccess().getCommaKeyword_6_0()); }
	','
	{ after(grammarAccess.getElkEdgeAccess().getCommaKeyword_6_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdge__Group_6__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdge__Group_6__1__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdge__Group_6__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeAccess().getTargetsAssignment_6_1()); }
	(rule__ElkEdge__TargetsAssignment_6_1)
	{ after(grammarAccess.getElkEdgeAccess().getTargetsAssignment_6_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkEdge__Group_7__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdge__Group_7__0__Impl
	rule__ElkEdge__Group_7__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdge__Group_7__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeAccess().getLeftCurlyBracketKeyword_7_0()); }
	'{'
	{ after(grammarAccess.getElkEdgeAccess().getLeftCurlyBracketKeyword_7_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdge__Group_7__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdge__Group_7__1__Impl
	rule__ElkEdge__Group_7__2
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdge__Group_7__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeAccess().getEdgeLayoutParserRuleCall_7_1()); }
	(ruleEdgeLayout)?
	{ after(grammarAccess.getElkEdgeAccess().getEdgeLayoutParserRuleCall_7_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdge__Group_7__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdge__Group_7__2__Impl
	rule__ElkEdge__Group_7__3
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdge__Group_7__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeAccess().getPropertiesAssignment_7_2()); }
	(rule__ElkEdge__PropertiesAssignment_7_2)*
	{ after(grammarAccess.getElkEdgeAccess().getPropertiesAssignment_7_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdge__Group_7__3
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdge__Group_7__3__Impl
	rule__ElkEdge__Group_7__4
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdge__Group_7__3__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeAccess().getLabelsAssignment_7_3()); }
	(rule__ElkEdge__LabelsAssignment_7_3)*
	{ after(grammarAccess.getElkEdgeAccess().getLabelsAssignment_7_3()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdge__Group_7__4
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdge__Group_7__4__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdge__Group_7__4__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeAccess().getRightCurlyBracketKeyword_7_4()); }
	'}'
	{ after(grammarAccess.getElkEdgeAccess().getRightCurlyBracketKeyword_7_4()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__EdgeLayout__Group__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__EdgeLayout__Group__0__Impl
	rule__EdgeLayout__Group__1
;
finally {
	restoreStackSize(stackSize);
}

rule__EdgeLayout__Group__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getEdgeLayoutAccess().getLayoutKeyword_0()); }
	'layout'
	{ after(grammarAccess.getEdgeLayoutAccess().getLayoutKeyword_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__EdgeLayout__Group__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__EdgeLayout__Group__1__Impl
	rule__EdgeLayout__Group__2
;
finally {
	restoreStackSize(stackSize);
}

rule__EdgeLayout__Group__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getEdgeLayoutAccess().getLeftSquareBracketKeyword_1()); }
	'['
	{ after(grammarAccess.getEdgeLayoutAccess().getLeftSquareBracketKeyword_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__EdgeLayout__Group__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__EdgeLayout__Group__2__Impl
	rule__EdgeLayout__Group__3
;
finally {
	restoreStackSize(stackSize);
}

rule__EdgeLayout__Group__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getEdgeLayoutAccess().getAlternatives_2()); }
	(rule__EdgeLayout__Alternatives_2)
	{ after(grammarAccess.getEdgeLayoutAccess().getAlternatives_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__EdgeLayout__Group__3
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__EdgeLayout__Group__3__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__EdgeLayout__Group__3__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getEdgeLayoutAccess().getRightSquareBracketKeyword_3()); }
	']'
	{ after(grammarAccess.getEdgeLayoutAccess().getRightSquareBracketKeyword_3()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkSingleEdgeSection__Group__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkSingleEdgeSection__Group__0__Impl
	rule__ElkSingleEdgeSection__Group__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__Group__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkSingleEdgeSectionAccess().getElkEdgeSectionAction_0()); }
	()
	{ after(grammarAccess.getElkSingleEdgeSectionAccess().getElkEdgeSectionAction_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__Group__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkSingleEdgeSection__Group__1__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__Group__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkSingleEdgeSectionAccess().getGroup_1()); }
	(rule__ElkSingleEdgeSection__Group_1__0)
	{ after(grammarAccess.getElkSingleEdgeSectionAccess().getGroup_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkSingleEdgeSection__Group_1__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkSingleEdgeSection__Group_1__0__Impl
	rule__ElkSingleEdgeSection__Group_1__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__Group_1__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0()); }
	(rule__ElkSingleEdgeSection__UnorderedGroup_1_0)
	{ after(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__Group_1__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkSingleEdgeSection__Group_1__1__Impl
	rule__ElkSingleEdgeSection__Group_1__2
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__Group_1__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkSingleEdgeSectionAccess().getGroup_1_1()); }
	(rule__ElkSingleEdgeSection__Group_1_1__0)?
	{ after(grammarAccess.getElkSingleEdgeSectionAccess().getGroup_1_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__Group_1__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkSingleEdgeSection__Group_1__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__Group_1__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkSingleEdgeSectionAccess().getPropertiesAssignment_1_2()); }
	(rule__ElkSingleEdgeSection__PropertiesAssignment_1_2)*
	{ after(grammarAccess.getElkSingleEdgeSectionAccess().getPropertiesAssignment_1_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkSingleEdgeSection__Group_1_0_0__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkSingleEdgeSection__Group_1_0_0__0__Impl
	rule__ElkSingleEdgeSection__Group_1_0_0__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__Group_1_0_0__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkSingleEdgeSectionAccess().getIncomingKeyword_1_0_0_0()); }
	'incoming'
	{ after(grammarAccess.getElkSingleEdgeSectionAccess().getIncomingKeyword_1_0_0_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__Group_1_0_0__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkSingleEdgeSection__Group_1_0_0__1__Impl
	rule__ElkSingleEdgeSection__Group_1_0_0__2
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__Group_1_0_0__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkSingleEdgeSectionAccess().getColonKeyword_1_0_0_1()); }
	':'
	{ after(grammarAccess.getElkSingleEdgeSectionAccess().getColonKeyword_1_0_0_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__Group_1_0_0__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkSingleEdgeSection__Group_1_0_0__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__Group_1_0_0__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkSingleEdgeSectionAccess().getIncomingShapeAssignment_1_0_0_2()); }
	(rule__ElkSingleEdgeSection__IncomingShapeAssignment_1_0_0_2)
	{ after(grammarAccess.getElkSingleEdgeSectionAccess().getIncomingShapeAssignment_1_0_0_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkSingleEdgeSection__Group_1_0_1__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkSingleEdgeSection__Group_1_0_1__0__Impl
	rule__ElkSingleEdgeSection__Group_1_0_1__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__Group_1_0_1__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkSingleEdgeSectionAccess().getOutgoingKeyword_1_0_1_0()); }
	'outgoing'
	{ after(grammarAccess.getElkSingleEdgeSectionAccess().getOutgoingKeyword_1_0_1_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__Group_1_0_1__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkSingleEdgeSection__Group_1_0_1__1__Impl
	rule__ElkSingleEdgeSection__Group_1_0_1__2
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__Group_1_0_1__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkSingleEdgeSectionAccess().getColonKeyword_1_0_1_1()); }
	':'
	{ after(grammarAccess.getElkSingleEdgeSectionAccess().getColonKeyword_1_0_1_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__Group_1_0_1__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkSingleEdgeSection__Group_1_0_1__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__Group_1_0_1__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkSingleEdgeSectionAccess().getOutgoingShapeAssignment_1_0_1_2()); }
	(rule__ElkSingleEdgeSection__OutgoingShapeAssignment_1_0_1_2)
	{ after(grammarAccess.getElkSingleEdgeSectionAccess().getOutgoingShapeAssignment_1_0_1_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkSingleEdgeSection__Group_1_0_2__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkSingleEdgeSection__Group_1_0_2__0__Impl
	rule__ElkSingleEdgeSection__Group_1_0_2__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__Group_1_0_2__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkSingleEdgeSectionAccess().getStartKeyword_1_0_2_0()); }
	'start'
	{ after(grammarAccess.getElkSingleEdgeSectionAccess().getStartKeyword_1_0_2_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__Group_1_0_2__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkSingleEdgeSection__Group_1_0_2__1__Impl
	rule__ElkSingleEdgeSection__Group_1_0_2__2
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__Group_1_0_2__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkSingleEdgeSectionAccess().getColonKeyword_1_0_2_1()); }
	':'
	{ after(grammarAccess.getElkSingleEdgeSectionAccess().getColonKeyword_1_0_2_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__Group_1_0_2__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkSingleEdgeSection__Group_1_0_2__2__Impl
	rule__ElkSingleEdgeSection__Group_1_0_2__3
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__Group_1_0_2__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkSingleEdgeSectionAccess().getStartXAssignment_1_0_2_2()); }
	(rule__ElkSingleEdgeSection__StartXAssignment_1_0_2_2)
	{ after(grammarAccess.getElkSingleEdgeSectionAccess().getStartXAssignment_1_0_2_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__Group_1_0_2__3
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkSingleEdgeSection__Group_1_0_2__3__Impl
	rule__ElkSingleEdgeSection__Group_1_0_2__4
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__Group_1_0_2__3__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkSingleEdgeSectionAccess().getCommaKeyword_1_0_2_3()); }
	','
	{ after(grammarAccess.getElkSingleEdgeSectionAccess().getCommaKeyword_1_0_2_3()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__Group_1_0_2__4
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkSingleEdgeSection__Group_1_0_2__4__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__Group_1_0_2__4__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkSingleEdgeSectionAccess().getStartYAssignment_1_0_2_4()); }
	(rule__ElkSingleEdgeSection__StartYAssignment_1_0_2_4)
	{ after(grammarAccess.getElkSingleEdgeSectionAccess().getStartYAssignment_1_0_2_4()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkSingleEdgeSection__Group_1_0_3__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkSingleEdgeSection__Group_1_0_3__0__Impl
	rule__ElkSingleEdgeSection__Group_1_0_3__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__Group_1_0_3__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkSingleEdgeSectionAccess().getEndKeyword_1_0_3_0()); }
	'end'
	{ after(grammarAccess.getElkSingleEdgeSectionAccess().getEndKeyword_1_0_3_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__Group_1_0_3__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkSingleEdgeSection__Group_1_0_3__1__Impl
	rule__ElkSingleEdgeSection__Group_1_0_3__2
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__Group_1_0_3__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkSingleEdgeSectionAccess().getColonKeyword_1_0_3_1()); }
	':'
	{ after(grammarAccess.getElkSingleEdgeSectionAccess().getColonKeyword_1_0_3_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__Group_1_0_3__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkSingleEdgeSection__Group_1_0_3__2__Impl
	rule__ElkSingleEdgeSection__Group_1_0_3__3
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__Group_1_0_3__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkSingleEdgeSectionAccess().getEndXAssignment_1_0_3_2()); }
	(rule__ElkSingleEdgeSection__EndXAssignment_1_0_3_2)
	{ after(grammarAccess.getElkSingleEdgeSectionAccess().getEndXAssignment_1_0_3_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__Group_1_0_3__3
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkSingleEdgeSection__Group_1_0_3__3__Impl
	rule__ElkSingleEdgeSection__Group_1_0_3__4
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__Group_1_0_3__3__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkSingleEdgeSectionAccess().getCommaKeyword_1_0_3_3()); }
	','
	{ after(grammarAccess.getElkSingleEdgeSectionAccess().getCommaKeyword_1_0_3_3()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__Group_1_0_3__4
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkSingleEdgeSection__Group_1_0_3__4__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__Group_1_0_3__4__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkSingleEdgeSectionAccess().getEndYAssignment_1_0_3_4()); }
	(rule__ElkSingleEdgeSection__EndYAssignment_1_0_3_4)
	{ after(grammarAccess.getElkSingleEdgeSectionAccess().getEndYAssignment_1_0_3_4()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkSingleEdgeSection__Group_1_1__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkSingleEdgeSection__Group_1_1__0__Impl
	rule__ElkSingleEdgeSection__Group_1_1__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__Group_1_1__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkSingleEdgeSectionAccess().getBendsKeyword_1_1_0()); }
	'bends'
	{ after(grammarAccess.getElkSingleEdgeSectionAccess().getBendsKeyword_1_1_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__Group_1_1__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkSingleEdgeSection__Group_1_1__1__Impl
	rule__ElkSingleEdgeSection__Group_1_1__2
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__Group_1_1__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkSingleEdgeSectionAccess().getColonKeyword_1_1_1()); }
	':'
	{ after(grammarAccess.getElkSingleEdgeSectionAccess().getColonKeyword_1_1_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__Group_1_1__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkSingleEdgeSection__Group_1_1__2__Impl
	rule__ElkSingleEdgeSection__Group_1_1__3
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__Group_1_1__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkSingleEdgeSectionAccess().getBendPointsAssignment_1_1_2()); }
	(rule__ElkSingleEdgeSection__BendPointsAssignment_1_1_2)
	{ after(grammarAccess.getElkSingleEdgeSectionAccess().getBendPointsAssignment_1_1_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__Group_1_1__3
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkSingleEdgeSection__Group_1_1__3__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__Group_1_1__3__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkSingleEdgeSectionAccess().getGroup_1_1_3()); }
	(rule__ElkSingleEdgeSection__Group_1_1_3__0)*
	{ after(grammarAccess.getElkSingleEdgeSectionAccess().getGroup_1_1_3()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkSingleEdgeSection__Group_1_1_3__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkSingleEdgeSection__Group_1_1_3__0__Impl
	rule__ElkSingleEdgeSection__Group_1_1_3__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__Group_1_1_3__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkSingleEdgeSectionAccess().getVerticalLineKeyword_1_1_3_0()); }
	'|'
	{ after(grammarAccess.getElkSingleEdgeSectionAccess().getVerticalLineKeyword_1_1_3_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__Group_1_1_3__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkSingleEdgeSection__Group_1_1_3__1__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__Group_1_1_3__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkSingleEdgeSectionAccess().getBendPointsAssignment_1_1_3_1()); }
	(rule__ElkSingleEdgeSection__BendPointsAssignment_1_1_3_1)
	{ after(grammarAccess.getElkSingleEdgeSectionAccess().getBendPointsAssignment_1_1_3_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkEdgeSection__Group__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeSection__Group__0__Impl
	rule__ElkEdgeSection__Group__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeSectionAccess().getSectionKeyword_0()); }
	'section'
	{ after(grammarAccess.getElkEdgeSectionAccess().getSectionKeyword_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeSection__Group__1__Impl
	rule__ElkEdgeSection__Group__2
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeSectionAccess().getIdentifierAssignment_1()); }
	(rule__ElkEdgeSection__IdentifierAssignment_1)
	{ after(grammarAccess.getElkEdgeSectionAccess().getIdentifierAssignment_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeSection__Group__2__Impl
	rule__ElkEdgeSection__Group__3
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeSectionAccess().getGroup_2()); }
	(rule__ElkEdgeSection__Group_2__0)?
	{ after(grammarAccess.getElkEdgeSectionAccess().getGroup_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group__3
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeSection__Group__3__Impl
	rule__ElkEdgeSection__Group__4
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group__3__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeSectionAccess().getLeftSquareBracketKeyword_3()); }
	'['
	{ after(grammarAccess.getElkEdgeSectionAccess().getLeftSquareBracketKeyword_3()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group__4
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeSection__Group__4__Impl
	rule__ElkEdgeSection__Group__5
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group__4__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeSectionAccess().getGroup_4()); }
	(rule__ElkEdgeSection__Group_4__0)
	{ after(grammarAccess.getElkEdgeSectionAccess().getGroup_4()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group__5
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeSection__Group__5__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group__5__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeSectionAccess().getRightSquareBracketKeyword_5()); }
	']'
	{ after(grammarAccess.getElkEdgeSectionAccess().getRightSquareBracketKeyword_5()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkEdgeSection__Group_2__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeSection__Group_2__0__Impl
	rule__ElkEdgeSection__Group_2__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group_2__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeSectionAccess().getHyphenMinusGreaterThanSignKeyword_2_0()); }
	'->'
	{ after(grammarAccess.getElkEdgeSectionAccess().getHyphenMinusGreaterThanSignKeyword_2_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group_2__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeSection__Group_2__1__Impl
	rule__ElkEdgeSection__Group_2__2
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group_2__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeSectionAccess().getOutgoingSectionsAssignment_2_1()); }
	(rule__ElkEdgeSection__OutgoingSectionsAssignment_2_1)
	{ after(grammarAccess.getElkEdgeSectionAccess().getOutgoingSectionsAssignment_2_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group_2__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeSection__Group_2__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group_2__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeSectionAccess().getGroup_2_2()); }
	(rule__ElkEdgeSection__Group_2_2__0)*
	{ after(grammarAccess.getElkEdgeSectionAccess().getGroup_2_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkEdgeSection__Group_2_2__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeSection__Group_2_2__0__Impl
	rule__ElkEdgeSection__Group_2_2__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group_2_2__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeSectionAccess().getCommaKeyword_2_2_0()); }
	','
	{ after(grammarAccess.getElkEdgeSectionAccess().getCommaKeyword_2_2_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group_2_2__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeSection__Group_2_2__1__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group_2_2__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeSectionAccess().getOutgoingSectionsAssignment_2_2_1()); }
	(rule__ElkEdgeSection__OutgoingSectionsAssignment_2_2_1)
	{ after(grammarAccess.getElkEdgeSectionAccess().getOutgoingSectionsAssignment_2_2_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkEdgeSection__Group_4__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeSection__Group_4__0__Impl
	rule__ElkEdgeSection__Group_4__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group_4__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0()); }
	(rule__ElkEdgeSection__UnorderedGroup_4_0)
	{ after(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group_4__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeSection__Group_4__1__Impl
	rule__ElkEdgeSection__Group_4__2
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group_4__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeSectionAccess().getGroup_4_1()); }
	(rule__ElkEdgeSection__Group_4_1__0)?
	{ after(grammarAccess.getElkEdgeSectionAccess().getGroup_4_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group_4__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeSection__Group_4__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group_4__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeSectionAccess().getPropertiesAssignment_4_2()); }
	(rule__ElkEdgeSection__PropertiesAssignment_4_2)*
	{ after(grammarAccess.getElkEdgeSectionAccess().getPropertiesAssignment_4_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkEdgeSection__Group_4_0_0__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeSection__Group_4_0_0__0__Impl
	rule__ElkEdgeSection__Group_4_0_0__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group_4_0_0__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeSectionAccess().getIncomingKeyword_4_0_0_0()); }
	'incoming'
	{ after(grammarAccess.getElkEdgeSectionAccess().getIncomingKeyword_4_0_0_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group_4_0_0__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeSection__Group_4_0_0__1__Impl
	rule__ElkEdgeSection__Group_4_0_0__2
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group_4_0_0__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeSectionAccess().getColonKeyword_4_0_0_1()); }
	':'
	{ after(grammarAccess.getElkEdgeSectionAccess().getColonKeyword_4_0_0_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group_4_0_0__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeSection__Group_4_0_0__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group_4_0_0__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeSectionAccess().getIncomingShapeAssignment_4_0_0_2()); }
	(rule__ElkEdgeSection__IncomingShapeAssignment_4_0_0_2)
	{ after(grammarAccess.getElkEdgeSectionAccess().getIncomingShapeAssignment_4_0_0_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkEdgeSection__Group_4_0_1__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeSection__Group_4_0_1__0__Impl
	rule__ElkEdgeSection__Group_4_0_1__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group_4_0_1__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeSectionAccess().getOutgoingKeyword_4_0_1_0()); }
	'outgoing'
	{ after(grammarAccess.getElkEdgeSectionAccess().getOutgoingKeyword_4_0_1_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group_4_0_1__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeSection__Group_4_0_1__1__Impl
	rule__ElkEdgeSection__Group_4_0_1__2
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group_4_0_1__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeSectionAccess().getColonKeyword_4_0_1_1()); }
	':'
	{ after(grammarAccess.getElkEdgeSectionAccess().getColonKeyword_4_0_1_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group_4_0_1__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeSection__Group_4_0_1__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group_4_0_1__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeSectionAccess().getOutgoingShapeAssignment_4_0_1_2()); }
	(rule__ElkEdgeSection__OutgoingShapeAssignment_4_0_1_2)
	{ after(grammarAccess.getElkEdgeSectionAccess().getOutgoingShapeAssignment_4_0_1_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkEdgeSection__Group_4_0_2__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeSection__Group_4_0_2__0__Impl
	rule__ElkEdgeSection__Group_4_0_2__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group_4_0_2__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeSectionAccess().getStartKeyword_4_0_2_0()); }
	'start'
	{ after(grammarAccess.getElkEdgeSectionAccess().getStartKeyword_4_0_2_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group_4_0_2__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeSection__Group_4_0_2__1__Impl
	rule__ElkEdgeSection__Group_4_0_2__2
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group_4_0_2__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeSectionAccess().getColonKeyword_4_0_2_1()); }
	':'
	{ after(grammarAccess.getElkEdgeSectionAccess().getColonKeyword_4_0_2_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group_4_0_2__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeSection__Group_4_0_2__2__Impl
	rule__ElkEdgeSection__Group_4_0_2__3
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group_4_0_2__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeSectionAccess().getStartXAssignment_4_0_2_2()); }
	(rule__ElkEdgeSection__StartXAssignment_4_0_2_2)
	{ after(grammarAccess.getElkEdgeSectionAccess().getStartXAssignment_4_0_2_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group_4_0_2__3
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeSection__Group_4_0_2__3__Impl
	rule__ElkEdgeSection__Group_4_0_2__4
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group_4_0_2__3__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeSectionAccess().getCommaKeyword_4_0_2_3()); }
	','
	{ after(grammarAccess.getElkEdgeSectionAccess().getCommaKeyword_4_0_2_3()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group_4_0_2__4
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeSection__Group_4_0_2__4__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group_4_0_2__4__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeSectionAccess().getStartYAssignment_4_0_2_4()); }
	(rule__ElkEdgeSection__StartYAssignment_4_0_2_4)
	{ after(grammarAccess.getElkEdgeSectionAccess().getStartYAssignment_4_0_2_4()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkEdgeSection__Group_4_0_3__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeSection__Group_4_0_3__0__Impl
	rule__ElkEdgeSection__Group_4_0_3__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group_4_0_3__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeSectionAccess().getEndKeyword_4_0_3_0()); }
	'end'
	{ after(grammarAccess.getElkEdgeSectionAccess().getEndKeyword_4_0_3_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group_4_0_3__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeSection__Group_4_0_3__1__Impl
	rule__ElkEdgeSection__Group_4_0_3__2
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group_4_0_3__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeSectionAccess().getColonKeyword_4_0_3_1()); }
	':'
	{ after(grammarAccess.getElkEdgeSectionAccess().getColonKeyword_4_0_3_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group_4_0_3__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeSection__Group_4_0_3__2__Impl
	rule__ElkEdgeSection__Group_4_0_3__3
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group_4_0_3__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeSectionAccess().getEndXAssignment_4_0_3_2()); }
	(rule__ElkEdgeSection__EndXAssignment_4_0_3_2)
	{ after(grammarAccess.getElkEdgeSectionAccess().getEndXAssignment_4_0_3_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group_4_0_3__3
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeSection__Group_4_0_3__3__Impl
	rule__ElkEdgeSection__Group_4_0_3__4
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group_4_0_3__3__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeSectionAccess().getCommaKeyword_4_0_3_3()); }
	','
	{ after(grammarAccess.getElkEdgeSectionAccess().getCommaKeyword_4_0_3_3()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group_4_0_3__4
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeSection__Group_4_0_3__4__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group_4_0_3__4__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeSectionAccess().getEndYAssignment_4_0_3_4()); }
	(rule__ElkEdgeSection__EndYAssignment_4_0_3_4)
	{ after(grammarAccess.getElkEdgeSectionAccess().getEndYAssignment_4_0_3_4()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkEdgeSection__Group_4_1__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeSection__Group_4_1__0__Impl
	rule__ElkEdgeSection__Group_4_1__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group_4_1__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeSectionAccess().getBendsKeyword_4_1_0()); }
	'bends'
	{ after(grammarAccess.getElkEdgeSectionAccess().getBendsKeyword_4_1_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group_4_1__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeSection__Group_4_1__1__Impl
	rule__ElkEdgeSection__Group_4_1__2
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group_4_1__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeSectionAccess().getColonKeyword_4_1_1()); }
	':'
	{ after(grammarAccess.getElkEdgeSectionAccess().getColonKeyword_4_1_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group_4_1__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeSection__Group_4_1__2__Impl
	rule__ElkEdgeSection__Group_4_1__3
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group_4_1__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeSectionAccess().getBendPointsAssignment_4_1_2()); }
	(rule__ElkEdgeSection__BendPointsAssignment_4_1_2)
	{ after(grammarAccess.getElkEdgeSectionAccess().getBendPointsAssignment_4_1_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group_4_1__3
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeSection__Group_4_1__3__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group_4_1__3__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeSectionAccess().getGroup_4_1_3()); }
	(rule__ElkEdgeSection__Group_4_1_3__0)*
	{ after(grammarAccess.getElkEdgeSectionAccess().getGroup_4_1_3()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkEdgeSection__Group_4_1_3__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeSection__Group_4_1_3__0__Impl
	rule__ElkEdgeSection__Group_4_1_3__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group_4_1_3__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeSectionAccess().getVerticalLineKeyword_4_1_3_0()); }
	'|'
	{ after(grammarAccess.getElkEdgeSectionAccess().getVerticalLineKeyword_4_1_3_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group_4_1_3__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeSection__Group_4_1_3__1__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__Group_4_1_3__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeSectionAccess().getBendPointsAssignment_4_1_3_1()); }
	(rule__ElkEdgeSection__BendPointsAssignment_4_1_3_1)
	{ after(grammarAccess.getElkEdgeSectionAccess().getBendPointsAssignment_4_1_3_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkBendPoint__Group__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkBendPoint__Group__0__Impl
	rule__ElkBendPoint__Group__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkBendPoint__Group__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkBendPointAccess().getXAssignment_0()); }
	(rule__ElkBendPoint__XAssignment_0)
	{ after(grammarAccess.getElkBendPointAccess().getXAssignment_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkBendPoint__Group__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkBendPoint__Group__1__Impl
	rule__ElkBendPoint__Group__2
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkBendPoint__Group__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkBendPointAccess().getCommaKeyword_1()); }
	','
	{ after(grammarAccess.getElkBendPointAccess().getCommaKeyword_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkBendPoint__Group__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkBendPoint__Group__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkBendPoint__Group__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkBendPointAccess().getYAssignment_2()); }
	(rule__ElkBendPoint__YAssignment_2)
	{ after(grammarAccess.getElkBendPointAccess().getYAssignment_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__QualifiedId__Group__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__QualifiedId__Group__0__Impl
	rule__QualifiedId__Group__1
;
finally {
	restoreStackSize(stackSize);
}

rule__QualifiedId__Group__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getQualifiedIdAccess().getIDTerminalRuleCall_0()); }
	RULE_ID
	{ after(grammarAccess.getQualifiedIdAccess().getIDTerminalRuleCall_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__QualifiedId__Group__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__QualifiedId__Group__1__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__QualifiedId__Group__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getQualifiedIdAccess().getGroup_1()); }
	(rule__QualifiedId__Group_1__0)*
	{ after(grammarAccess.getQualifiedIdAccess().getGroup_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__QualifiedId__Group_1__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__QualifiedId__Group_1__0__Impl
	rule__QualifiedId__Group_1__1
;
finally {
	restoreStackSize(stackSize);
}

rule__QualifiedId__Group_1__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getQualifiedIdAccess().getFullStopKeyword_1_0()); }
	'.'
	{ after(grammarAccess.getQualifiedIdAccess().getFullStopKeyword_1_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__QualifiedId__Group_1__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__QualifiedId__Group_1__1__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__QualifiedId__Group_1__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getQualifiedIdAccess().getIDTerminalRuleCall_1_1()); }
	RULE_ID
	{ after(grammarAccess.getQualifiedIdAccess().getIDTerminalRuleCall_1_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__Property__Group__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Property__Group__0__Impl
	rule__Property__Group__1
;
finally {
	restoreStackSize(stackSize);
}

rule__Property__Group__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getPropertyAccess().getKeyAssignment_0()); }
	(rule__Property__KeyAssignment_0)
	{ after(grammarAccess.getPropertyAccess().getKeyAssignment_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Property__Group__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Property__Group__1__Impl
	rule__Property__Group__2
;
finally {
	restoreStackSize(stackSize);
}

rule__Property__Group__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getPropertyAccess().getColonKeyword_1()); }
	':'
	{ after(grammarAccess.getPropertyAccess().getColonKeyword_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__Property__Group__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__Property__Group__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__Property__Group__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getPropertyAccess().getAlternatives_2()); }
	(rule__Property__Alternatives_2)
	{ after(grammarAccess.getPropertyAccess().getAlternatives_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__PropertyKey__Group__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__PropertyKey__Group__0__Impl
	rule__PropertyKey__Group__1
;
finally {
	restoreStackSize(stackSize);
}

rule__PropertyKey__Group__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getPropertyKeyAccess().getIDTerminalRuleCall_0()); }
	RULE_ID
	{ after(grammarAccess.getPropertyKeyAccess().getIDTerminalRuleCall_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__PropertyKey__Group__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__PropertyKey__Group__1__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__PropertyKey__Group__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getPropertyKeyAccess().getGroup_1()); }
	(rule__PropertyKey__Group_1__0)*
	{ after(grammarAccess.getPropertyKeyAccess().getGroup_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__PropertyKey__Group_1__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__PropertyKey__Group_1__0__Impl
	rule__PropertyKey__Group_1__1
;
finally {
	restoreStackSize(stackSize);
}

rule__PropertyKey__Group_1__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getPropertyKeyAccess().getFullStopKeyword_1_0()); }
	'.'
	{ after(grammarAccess.getPropertyKeyAccess().getFullStopKeyword_1_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__PropertyKey__Group_1__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__PropertyKey__Group_1__1__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__PropertyKey__Group_1__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getPropertyKeyAccess().getIDTerminalRuleCall_1_1()); }
	RULE_ID
	{ after(grammarAccess.getPropertyKeyAccess().getIDTerminalRuleCall_1_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ShapeLayout__UnorderedGroup_2
	@init {
		int stackSize = keepStackSize();
		getUnorderedGroupHelper().enter(grammarAccess.getShapeLayoutAccess().getUnorderedGroup_2());
	}
:
	rule__ShapeLayout__UnorderedGroup_2__0
	?
;
finally {
	getUnorderedGroupHelper().leave(grammarAccess.getShapeLayoutAccess().getUnorderedGroup_2());
	restoreStackSize(stackSize);
}

rule__ShapeLayout__UnorderedGroup_2__Impl
	@init {
		int stackSize = keepStackSize();
		boolean selected = false;
	}
:
		(
		( 
			{getUnorderedGroupHelper().canSelect(grammarAccess.getShapeLayoutAccess().getUnorderedGroup_2(), 0)}?=>(
				{
					getUnorderedGroupHelper().select(grammarAccess.getShapeLayoutAccess().getUnorderedGroup_2(), 0);
				}
				{
					selected = true;
				}
				(
					{ before(grammarAccess.getShapeLayoutAccess().getGroup_2_0()); }
					(rule__ShapeLayout__Group_2_0__0)
					{ after(grammarAccess.getShapeLayoutAccess().getGroup_2_0()); }
				)
			)
		)|
		( 
			{getUnorderedGroupHelper().canSelect(grammarAccess.getShapeLayoutAccess().getUnorderedGroup_2(), 1)}?=>(
				{
					getUnorderedGroupHelper().select(grammarAccess.getShapeLayoutAccess().getUnorderedGroup_2(), 1);
				}
				{
					selected = true;
				}
				(
					{ before(grammarAccess.getShapeLayoutAccess().getGroup_2_1()); }
					(rule__ShapeLayout__Group_2_1__0)
					{ after(grammarAccess.getShapeLayoutAccess().getGroup_2_1()); }
				)
			)
		)
		)
;
finally {
	if (selected)
		getUnorderedGroupHelper().returnFromSelection(grammarAccess.getShapeLayoutAccess().getUnorderedGroup_2());
	restoreStackSize(stackSize);
}

rule__ShapeLayout__UnorderedGroup_2__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ShapeLayout__UnorderedGroup_2__Impl
	rule__ShapeLayout__UnorderedGroup_2__1?
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeLayout__UnorderedGroup_2__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ShapeLayout__UnorderedGroup_2__Impl
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkSingleEdgeSection__UnorderedGroup_1_0
	@init {
		int stackSize = keepStackSize();
		getUnorderedGroupHelper().enter(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0());
	}
:
	rule__ElkSingleEdgeSection__UnorderedGroup_1_0__0
	?
;
finally {
	getUnorderedGroupHelper().leave(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0());
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__UnorderedGroup_1_0__Impl
	@init {
		int stackSize = keepStackSize();
		boolean selected = false;
	}
:
		(
		( 
			{getUnorderedGroupHelper().canSelect(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0(), 0)}?=>(
				{
					getUnorderedGroupHelper().select(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0(), 0);
				}
				{
					selected = true;
				}
				(
					{ before(grammarAccess.getElkSingleEdgeSectionAccess().getGroup_1_0_0()); }
					(rule__ElkSingleEdgeSection__Group_1_0_0__0)
					{ after(grammarAccess.getElkSingleEdgeSectionAccess().getGroup_1_0_0()); }
				)
			)
		)|
		( 
			{getUnorderedGroupHelper().canSelect(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0(), 1)}?=>(
				{
					getUnorderedGroupHelper().select(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0(), 1);
				}
				{
					selected = true;
				}
				(
					{ before(grammarAccess.getElkSingleEdgeSectionAccess().getGroup_1_0_1()); }
					(rule__ElkSingleEdgeSection__Group_1_0_1__0)
					{ after(grammarAccess.getElkSingleEdgeSectionAccess().getGroup_1_0_1()); }
				)
			)
		)|
		( 
			{getUnorderedGroupHelper().canSelect(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0(), 2)}?=>(
				{
					getUnorderedGroupHelper().select(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0(), 2);
				}
				{
					selected = true;
				}
				(
					{ before(grammarAccess.getElkSingleEdgeSectionAccess().getGroup_1_0_2()); }
					(rule__ElkSingleEdgeSection__Group_1_0_2__0)
					{ after(grammarAccess.getElkSingleEdgeSectionAccess().getGroup_1_0_2()); }
				)
			)
		)|
		( 
			{getUnorderedGroupHelper().canSelect(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0(), 3)}?=>(
				{
					getUnorderedGroupHelper().select(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0(), 3);
				}
				{
					selected = true;
				}
				(
					{ before(grammarAccess.getElkSingleEdgeSectionAccess().getGroup_1_0_3()); }
					(rule__ElkSingleEdgeSection__Group_1_0_3__0)
					{ after(grammarAccess.getElkSingleEdgeSectionAccess().getGroup_1_0_3()); }
				)
			)
		)
		)
;
finally {
	if (selected)
		getUnorderedGroupHelper().returnFromSelection(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0());
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__UnorderedGroup_1_0__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkSingleEdgeSection__UnorderedGroup_1_0__Impl
	rule__ElkSingleEdgeSection__UnorderedGroup_1_0__1?
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__UnorderedGroup_1_0__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkSingleEdgeSection__UnorderedGroup_1_0__Impl
	rule__ElkSingleEdgeSection__UnorderedGroup_1_0__2?
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__UnorderedGroup_1_0__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkSingleEdgeSection__UnorderedGroup_1_0__Impl
	rule__ElkSingleEdgeSection__UnorderedGroup_1_0__3?
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__UnorderedGroup_1_0__3
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkSingleEdgeSection__UnorderedGroup_1_0__Impl
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkEdgeSection__UnorderedGroup_4_0
	@init {
		int stackSize = keepStackSize();
		getUnorderedGroupHelper().enter(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0());
	}
:
	rule__ElkEdgeSection__UnorderedGroup_4_0__0
	?
;
finally {
	getUnorderedGroupHelper().leave(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0());
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__UnorderedGroup_4_0__Impl
	@init {
		int stackSize = keepStackSize();
		boolean selected = false;
	}
:
		(
		( 
			{getUnorderedGroupHelper().canSelect(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0(), 0)}?=>(
				{
					getUnorderedGroupHelper().select(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0(), 0);
				}
				{
					selected = true;
				}
				(
					{ before(grammarAccess.getElkEdgeSectionAccess().getGroup_4_0_0()); }
					(rule__ElkEdgeSection__Group_4_0_0__0)
					{ after(grammarAccess.getElkEdgeSectionAccess().getGroup_4_0_0()); }
				)
			)
		)|
		( 
			{getUnorderedGroupHelper().canSelect(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0(), 1)}?=>(
				{
					getUnorderedGroupHelper().select(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0(), 1);
				}
				{
					selected = true;
				}
				(
					{ before(grammarAccess.getElkEdgeSectionAccess().getGroup_4_0_1()); }
					(rule__ElkEdgeSection__Group_4_0_1__0)
					{ after(grammarAccess.getElkEdgeSectionAccess().getGroup_4_0_1()); }
				)
			)
		)|
		( 
			{getUnorderedGroupHelper().canSelect(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0(), 2)}?=>(
				{
					getUnorderedGroupHelper().select(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0(), 2);
				}
				{
					selected = true;
				}
				(
					{ before(grammarAccess.getElkEdgeSectionAccess().getGroup_4_0_2()); }
					(rule__ElkEdgeSection__Group_4_0_2__0)
					{ after(grammarAccess.getElkEdgeSectionAccess().getGroup_4_0_2()); }
				)
			)
		)|
		( 
			{getUnorderedGroupHelper().canSelect(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0(), 3)}?=>(
				{
					getUnorderedGroupHelper().select(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0(), 3);
				}
				{
					selected = true;
				}
				(
					{ before(grammarAccess.getElkEdgeSectionAccess().getGroup_4_0_3()); }
					(rule__ElkEdgeSection__Group_4_0_3__0)
					{ after(grammarAccess.getElkEdgeSectionAccess().getGroup_4_0_3()); }
				)
			)
		)
		)
;
finally {
	if (selected)
		getUnorderedGroupHelper().returnFromSelection(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0());
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__UnorderedGroup_4_0__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeSection__UnorderedGroup_4_0__Impl
	rule__ElkEdgeSection__UnorderedGroup_4_0__1?
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__UnorderedGroup_4_0__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeSection__UnorderedGroup_4_0__Impl
	rule__ElkEdgeSection__UnorderedGroup_4_0__2?
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__UnorderedGroup_4_0__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeSection__UnorderedGroup_4_0__Impl
	rule__ElkEdgeSection__UnorderedGroup_4_0__3?
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__UnorderedGroup_4_0__3
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeSection__UnorderedGroup_4_0__Impl
;
finally {
	restoreStackSize(stackSize);
}


rule__RootNode__IdentifierAssignment_1_1
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getRootNodeAccess().getIdentifierIDTerminalRuleCall_1_1_0()); }
		RULE_ID
		{ after(grammarAccess.getRootNodeAccess().getIdentifierIDTerminalRuleCall_1_1_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__RootNode__PropertiesAssignment_3
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getRootNodeAccess().getPropertiesPropertyParserRuleCall_3_0()); }
		ruleProperty
		{ after(grammarAccess.getRootNodeAccess().getPropertiesPropertyParserRuleCall_3_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__RootNode__LabelsAssignment_4_0
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getRootNodeAccess().getLabelsElkLabelParserRuleCall_4_0_0()); }
		ruleElkLabel
		{ after(grammarAccess.getRootNodeAccess().getLabelsElkLabelParserRuleCall_4_0_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__RootNode__PortsAssignment_4_1
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getRootNodeAccess().getPortsElkPortParserRuleCall_4_1_0()); }
		ruleElkPort
		{ after(grammarAccess.getRootNodeAccess().getPortsElkPortParserRuleCall_4_1_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__RootNode__ChildrenAssignment_4_2
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getRootNodeAccess().getChildrenElkNodeParserRuleCall_4_2_0()); }
		ruleElkNode
		{ after(grammarAccess.getRootNodeAccess().getChildrenElkNodeParserRuleCall_4_2_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__RootNode__ContainedEdgesAssignment_4_3
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getRootNodeAccess().getContainedEdgesElkEdgeParserRuleCall_4_3_0()); }
		ruleElkEdge
		{ after(grammarAccess.getRootNodeAccess().getContainedEdgesElkEdgeParserRuleCall_4_3_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNode__IdentifierAssignment_1
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkNodeAccess().getIdentifierIDTerminalRuleCall_1_0()); }
		RULE_ID
		{ after(grammarAccess.getElkNodeAccess().getIdentifierIDTerminalRuleCall_1_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNode__PropertiesAssignment_2_2
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkNodeAccess().getPropertiesPropertyParserRuleCall_2_2_0()); }
		ruleProperty
		{ after(grammarAccess.getElkNodeAccess().getPropertiesPropertyParserRuleCall_2_2_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNode__LabelsAssignment_2_3_0
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkNodeAccess().getLabelsElkLabelParserRuleCall_2_3_0_0()); }
		ruleElkLabel
		{ after(grammarAccess.getElkNodeAccess().getLabelsElkLabelParserRuleCall_2_3_0_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNode__PortsAssignment_2_3_1
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkNodeAccess().getPortsElkPortParserRuleCall_2_3_1_0()); }
		ruleElkPort
		{ after(grammarAccess.getElkNodeAccess().getPortsElkPortParserRuleCall_2_3_1_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNode__ChildrenAssignment_2_3_2
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkNodeAccess().getChildrenElkNodeParserRuleCall_2_3_2_0()); }
		ruleElkNode
		{ after(grammarAccess.getElkNodeAccess().getChildrenElkNodeParserRuleCall_2_3_2_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNode__ContainedEdgesAssignment_2_3_3
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkNodeAccess().getContainedEdgesElkEdgeParserRuleCall_2_3_3_0()); }
		ruleElkEdge
		{ after(grammarAccess.getElkNodeAccess().getContainedEdgesElkEdgeParserRuleCall_2_3_3_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkLabel__IdentifierAssignment_1_0
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkLabelAccess().getIdentifierIDTerminalRuleCall_1_0_0()); }
		RULE_ID
		{ after(grammarAccess.getElkLabelAccess().getIdentifierIDTerminalRuleCall_1_0_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkLabel__TextAssignment_2
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkLabelAccess().getTextSTRINGTerminalRuleCall_2_0()); }
		RULE_STRING
		{ after(grammarAccess.getElkLabelAccess().getTextSTRINGTerminalRuleCall_2_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkLabel__PropertiesAssignment_3_2
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkLabelAccess().getPropertiesPropertyParserRuleCall_3_2_0()); }
		ruleProperty
		{ after(grammarAccess.getElkLabelAccess().getPropertiesPropertyParserRuleCall_3_2_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkLabel__LabelsAssignment_3_3
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkLabelAccess().getLabelsElkLabelParserRuleCall_3_3_0()); }
		ruleElkLabel
		{ after(grammarAccess.getElkLabelAccess().getLabelsElkLabelParserRuleCall_3_3_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkPort__IdentifierAssignment_1
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkPortAccess().getIdentifierIDTerminalRuleCall_1_0()); }
		RULE_ID
		{ after(grammarAccess.getElkPortAccess().getIdentifierIDTerminalRuleCall_1_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkPort__PropertiesAssignment_2_2
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkPortAccess().getPropertiesPropertyParserRuleCall_2_2_0()); }
		ruleProperty
		{ after(grammarAccess.getElkPortAccess().getPropertiesPropertyParserRuleCall_2_2_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkPort__LabelsAssignment_2_3
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkPortAccess().getLabelsElkLabelParserRuleCall_2_3_0()); }
		ruleElkLabel
		{ after(grammarAccess.getElkPortAccess().getLabelsElkLabelParserRuleCall_2_3_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeLayout__XAssignment_2_0_2
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getShapeLayoutAccess().getXNumberParserRuleCall_2_0_2_0()); }
		ruleNumber
		{ after(grammarAccess.getShapeLayoutAccess().getXNumberParserRuleCall_2_0_2_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeLayout__YAssignment_2_0_4
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getShapeLayoutAccess().getYNumberParserRuleCall_2_0_4_0()); }
		ruleNumber
		{ after(grammarAccess.getShapeLayoutAccess().getYNumberParserRuleCall_2_0_4_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeLayout__WidthAssignment_2_1_2
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getShapeLayoutAccess().getWidthNumberParserRuleCall_2_1_2_0()); }
		ruleNumber
		{ after(grammarAccess.getShapeLayoutAccess().getWidthNumberParserRuleCall_2_1_2_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeLayout__HeightAssignment_2_1_4
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getShapeLayoutAccess().getHeightNumberParserRuleCall_2_1_4_0()); }
		ruleNumber
		{ after(grammarAccess.getShapeLayoutAccess().getHeightNumberParserRuleCall_2_1_4_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdge__IdentifierAssignment_1_0
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkEdgeAccess().getIdentifierIDTerminalRuleCall_1_0_0()); }
		RULE_ID
		{ after(grammarAccess.getElkEdgeAccess().getIdentifierIDTerminalRuleCall_1_0_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdge__SourcesAssignment_2
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkEdgeAccess().getSourcesElkConnectableShapeCrossReference_2_0()); }
		(
			{ before(grammarAccess.getElkEdgeAccess().getSourcesElkConnectableShapeQualifiedIdParserRuleCall_2_0_1()); }
			ruleQualifiedId
			{ after(grammarAccess.getElkEdgeAccess().getSourcesElkConnectableShapeQualifiedIdParserRuleCall_2_0_1()); }
		)
		{ after(grammarAccess.getElkEdgeAccess().getSourcesElkConnectableShapeCrossReference_2_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdge__SourcesAssignment_3_1
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkEdgeAccess().getSourcesElkConnectableShapeCrossReference_3_1_0()); }
		(
			{ before(grammarAccess.getElkEdgeAccess().getSourcesElkConnectableShapeQualifiedIdParserRuleCall_3_1_0_1()); }
			ruleQualifiedId
			{ after(grammarAccess.getElkEdgeAccess().getSourcesElkConnectableShapeQualifiedIdParserRuleCall_3_1_0_1()); }
		)
		{ after(grammarAccess.getElkEdgeAccess().getSourcesElkConnectableShapeCrossReference_3_1_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdge__TargetsAssignment_5
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkEdgeAccess().getTargetsElkConnectableShapeCrossReference_5_0()); }
		(
			{ before(grammarAccess.getElkEdgeAccess().getTargetsElkConnectableShapeQualifiedIdParserRuleCall_5_0_1()); }
			ruleQualifiedId
			{ after(grammarAccess.getElkEdgeAccess().getTargetsElkConnectableShapeQualifiedIdParserRuleCall_5_0_1()); }
		)
		{ after(grammarAccess.getElkEdgeAccess().getTargetsElkConnectableShapeCrossReference_5_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdge__TargetsAssignment_6_1
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkEdgeAccess().getTargetsElkConnectableShapeCrossReference_6_1_0()); }
		(
			{ before(grammarAccess.getElkEdgeAccess().getTargetsElkConnectableShapeQualifiedIdParserRuleCall_6_1_0_1()); }
			ruleQualifiedId
			{ after(grammarAccess.getElkEdgeAccess().getTargetsElkConnectableShapeQualifiedIdParserRuleCall_6_1_0_1()); }
		)
		{ after(grammarAccess.getElkEdgeAccess().getTargetsElkConnectableShapeCrossReference_6_1_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdge__PropertiesAssignment_7_2
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkEdgeAccess().getPropertiesPropertyParserRuleCall_7_2_0()); }
		ruleProperty
		{ after(grammarAccess.getElkEdgeAccess().getPropertiesPropertyParserRuleCall_7_2_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdge__LabelsAssignment_7_3
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkEdgeAccess().getLabelsElkLabelParserRuleCall_7_3_0()); }
		ruleElkLabel
		{ after(grammarAccess.getElkEdgeAccess().getLabelsElkLabelParserRuleCall_7_3_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__EdgeLayout__SectionsAssignment_2_0
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getEdgeLayoutAccess().getSectionsElkSingleEdgeSectionParserRuleCall_2_0_0()); }
		ruleElkSingleEdgeSection
		{ after(grammarAccess.getEdgeLayoutAccess().getSectionsElkSingleEdgeSectionParserRuleCall_2_0_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__EdgeLayout__SectionsAssignment_2_1
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getEdgeLayoutAccess().getSectionsElkEdgeSectionParserRuleCall_2_1_0()); }
		ruleElkEdgeSection
		{ after(grammarAccess.getEdgeLayoutAccess().getSectionsElkEdgeSectionParserRuleCall_2_1_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__IncomingShapeAssignment_1_0_0_2
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkSingleEdgeSectionAccess().getIncomingShapeElkConnectableShapeCrossReference_1_0_0_2_0()); }
		(
			{ before(grammarAccess.getElkSingleEdgeSectionAccess().getIncomingShapeElkConnectableShapeQualifiedIdParserRuleCall_1_0_0_2_0_1()); }
			ruleQualifiedId
			{ after(grammarAccess.getElkSingleEdgeSectionAccess().getIncomingShapeElkConnectableShapeQualifiedIdParserRuleCall_1_0_0_2_0_1()); }
		)
		{ after(grammarAccess.getElkSingleEdgeSectionAccess().getIncomingShapeElkConnectableShapeCrossReference_1_0_0_2_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__OutgoingShapeAssignment_1_0_1_2
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkSingleEdgeSectionAccess().getOutgoingShapeElkConnectableShapeCrossReference_1_0_1_2_0()); }
		(
			{ before(grammarAccess.getElkSingleEdgeSectionAccess().getOutgoingShapeElkConnectableShapeQualifiedIdParserRuleCall_1_0_1_2_0_1()); }
			ruleQualifiedId
			{ after(grammarAccess.getElkSingleEdgeSectionAccess().getOutgoingShapeElkConnectableShapeQualifiedIdParserRuleCall_1_0_1_2_0_1()); }
		)
		{ after(grammarAccess.getElkSingleEdgeSectionAccess().getOutgoingShapeElkConnectableShapeCrossReference_1_0_1_2_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__StartXAssignment_1_0_2_2
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkSingleEdgeSectionAccess().getStartXNumberParserRuleCall_1_0_2_2_0()); }
		ruleNumber
		{ after(grammarAccess.getElkSingleEdgeSectionAccess().getStartXNumberParserRuleCall_1_0_2_2_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__StartYAssignment_1_0_2_4
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkSingleEdgeSectionAccess().getStartYNumberParserRuleCall_1_0_2_4_0()); }
		ruleNumber
		{ after(grammarAccess.getElkSingleEdgeSectionAccess().getStartYNumberParserRuleCall_1_0_2_4_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__EndXAssignment_1_0_3_2
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkSingleEdgeSectionAccess().getEndXNumberParserRuleCall_1_0_3_2_0()); }
		ruleNumber
		{ after(grammarAccess.getElkSingleEdgeSectionAccess().getEndXNumberParserRuleCall_1_0_3_2_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__EndYAssignment_1_0_3_4
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkSingleEdgeSectionAccess().getEndYNumberParserRuleCall_1_0_3_4_0()); }
		ruleNumber
		{ after(grammarAccess.getElkSingleEdgeSectionAccess().getEndYNumberParserRuleCall_1_0_3_4_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__BendPointsAssignment_1_1_2
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkSingleEdgeSectionAccess().getBendPointsElkBendPointParserRuleCall_1_1_2_0()); }
		ruleElkBendPoint
		{ after(grammarAccess.getElkSingleEdgeSectionAccess().getBendPointsElkBendPointParserRuleCall_1_1_2_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__BendPointsAssignment_1_1_3_1
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkSingleEdgeSectionAccess().getBendPointsElkBendPointParserRuleCall_1_1_3_1_0()); }
		ruleElkBendPoint
		{ after(grammarAccess.getElkSingleEdgeSectionAccess().getBendPointsElkBendPointParserRuleCall_1_1_3_1_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkSingleEdgeSection__PropertiesAssignment_1_2
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkSingleEdgeSectionAccess().getPropertiesPropertyParserRuleCall_1_2_0()); }
		ruleProperty
		{ after(grammarAccess.getElkSingleEdgeSectionAccess().getPropertiesPropertyParserRuleCall_1_2_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__IdentifierAssignment_1
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkEdgeSectionAccess().getIdentifierIDTerminalRuleCall_1_0()); }
		RULE_ID
		{ after(grammarAccess.getElkEdgeSectionAccess().getIdentifierIDTerminalRuleCall_1_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__OutgoingSectionsAssignment_2_1
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkEdgeSectionAccess().getOutgoingSectionsElkEdgeSectionCrossReference_2_1_0()); }
		(
			{ before(grammarAccess.getElkEdgeSectionAccess().getOutgoingSectionsElkEdgeSectionIDTerminalRuleCall_2_1_0_1()); }
			RULE_ID
			{ after(grammarAccess.getElkEdgeSectionAccess().getOutgoingSectionsElkEdgeSectionIDTerminalRuleCall_2_1_0_1()); }
		)
		{ after(grammarAccess.getElkEdgeSectionAccess().getOutgoingSectionsElkEdgeSectionCrossReference_2_1_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__OutgoingSectionsAssignment_2_2_1
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkEdgeSectionAccess().getOutgoingSectionsElkEdgeSectionCrossReference_2_2_1_0()); }
		(
			{ before(grammarAccess.getElkEdgeSectionAccess().getOutgoingSectionsElkEdgeSectionIDTerminalRuleCall_2_2_1_0_1()); }
			RULE_ID
			{ after(grammarAccess.getElkEdgeSectionAccess().getOutgoingSectionsElkEdgeSectionIDTerminalRuleCall_2_2_1_0_1()); }
		)
		{ after(grammarAccess.getElkEdgeSectionAccess().getOutgoingSectionsElkEdgeSectionCrossReference_2_2_1_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__IncomingShapeAssignment_4_0_0_2
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkEdgeSectionAccess().getIncomingShapeElkConnectableShapeCrossReference_4_0_0_2_0()); }
		(
			{ before(grammarAccess.getElkEdgeSectionAccess().getIncomingShapeElkConnectableShapeQualifiedIdParserRuleCall_4_0_0_2_0_1()); }
			ruleQualifiedId
			{ after(grammarAccess.getElkEdgeSectionAccess().getIncomingShapeElkConnectableShapeQualifiedIdParserRuleCall_4_0_0_2_0_1()); }
		)
		{ after(grammarAccess.getElkEdgeSectionAccess().getIncomingShapeElkConnectableShapeCrossReference_4_0_0_2_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__OutgoingShapeAssignment_4_0_1_2
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkEdgeSectionAccess().getOutgoingShapeElkConnectableShapeCrossReference_4_0_1_2_0()); }
		(
			{ before(grammarAccess.getElkEdgeSectionAccess().getOutgoingShapeElkConnectableShapeQualifiedIdParserRuleCall_4_0_1_2_0_1()); }
			ruleQualifiedId
			{ after(grammarAccess.getElkEdgeSectionAccess().getOutgoingShapeElkConnectableShapeQualifiedIdParserRuleCall_4_0_1_2_0_1()); }
		)
		{ after(grammarAccess.getElkEdgeSectionAccess().getOutgoingShapeElkConnectableShapeCrossReference_4_0_1_2_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__StartXAssignment_4_0_2_2
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkEdgeSectionAccess().getStartXNumberParserRuleCall_4_0_2_2_0()); }
		ruleNumber
		{ after(grammarAccess.getElkEdgeSectionAccess().getStartXNumberParserRuleCall_4_0_2_2_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__StartYAssignment_4_0_2_4
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkEdgeSectionAccess().getStartYNumberParserRuleCall_4_0_2_4_0()); }
		ruleNumber
		{ after(grammarAccess.getElkEdgeSectionAccess().getStartYNumberParserRuleCall_4_0_2_4_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__EndXAssignment_4_0_3_2
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkEdgeSectionAccess().getEndXNumberParserRuleCall_4_0_3_2_0()); }
		ruleNumber
		{ after(grammarAccess.getElkEdgeSectionAccess().getEndXNumberParserRuleCall_4_0_3_2_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__EndYAssignment_4_0_3_4
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkEdgeSectionAccess().getEndYNumberParserRuleCall_4_0_3_4_0()); }
		ruleNumber
		{ after(grammarAccess.getElkEdgeSectionAccess().getEndYNumberParserRuleCall_4_0_3_4_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__BendPointsAssignment_4_1_2
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkEdgeSectionAccess().getBendPointsElkBendPointParserRuleCall_4_1_2_0()); }
		ruleElkBendPoint
		{ after(grammarAccess.getElkEdgeSectionAccess().getBendPointsElkBendPointParserRuleCall_4_1_2_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__BendPointsAssignment_4_1_3_1
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkEdgeSectionAccess().getBendPointsElkBendPointParserRuleCall_4_1_3_1_0()); }
		ruleElkBendPoint
		{ after(grammarAccess.getElkEdgeSectionAccess().getBendPointsElkBendPointParserRuleCall_4_1_3_1_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSection__PropertiesAssignment_4_2
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkEdgeSectionAccess().getPropertiesPropertyParserRuleCall_4_2_0()); }
		ruleProperty
		{ after(grammarAccess.getElkEdgeSectionAccess().getPropertiesPropertyParserRuleCall_4_2_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkBendPoint__XAssignment_0
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkBendPointAccess().getXNumberParserRuleCall_0_0()); }
		ruleNumber
		{ after(grammarAccess.getElkBendPointAccess().getXNumberParserRuleCall_0_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkBendPoint__YAssignment_2
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkBendPointAccess().getYNumberParserRuleCall_2_0()); }
		ruleNumber
		{ after(grammarAccess.getElkBendPointAccess().getYNumberParserRuleCall_2_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__Property__KeyAssignment_0
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getPropertyAccess().getKeyPropertyKeyParserRuleCall_0_0()); }
		rulePropertyKey
		{ after(grammarAccess.getPropertyAccess().getKeyPropertyKeyParserRuleCall_0_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__Property__ValueAssignment_2_0
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getPropertyAccess().getValueStringValueParserRuleCall_2_0_0()); }
		ruleStringValue
		{ after(grammarAccess.getPropertyAccess().getValueStringValueParserRuleCall_2_0_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__Property__ValueAssignment_2_1
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getPropertyAccess().getValueQualifiedIdValueParserRuleCall_2_1_0()); }
		ruleQualifiedIdValue
		{ after(grammarAccess.getPropertyAccess().getValueQualifiedIdValueParserRuleCall_2_1_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__Property__ValueAssignment_2_2
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getPropertyAccess().getValueNumberValueParserRuleCall_2_2_0()); }
		ruleNumberValue
		{ after(grammarAccess.getPropertyAccess().getValueNumberValueParserRuleCall_2_2_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__Property__ValueAssignment_2_3
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getPropertyAccess().getValueBooleanValueParserRuleCall_2_3_0()); }
		ruleBooleanValue
		{ after(grammarAccess.getPropertyAccess().getValueBooleanValueParserRuleCall_2_3_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

RULE_SIGNED_INT : ('+'|'-')? RULE_INT;

RULE_FLOAT : ('+'|'-')? (RULE_INT '.' RULE_INT|RULE_INT ('.' RULE_INT)? ('e'|'E') ('+'|'-')? RULE_INT);

RULE_ID : '^'? ('a'..'z'|'A'..'Z'|'_') ('a'..'z'|'A'..'Z'|'_'|'0'..'9')*;

fragment RULE_INT : ('0'..'9')+;

RULE_STRING : ('"' ('\\' .|~(('\\'|'"')))* '"'|'\'' ('\\' .|~(('\\'|'\'')))* '\'');

RULE_ML_COMMENT : '/*' ( options {greedy=false;} : . )*'*/';

RULE_SL_COMMENT : '//' ~(('\n'|'\r'))* ('\r'? '\n')?;

RULE_WS : (' '|'\t'|'\r'|'\n')+;

RULE_ANY_OTHER : .;
