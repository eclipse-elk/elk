/*******************************************************************************
 * Copyright (c) 2020 Kiel University and others.
 *
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * This Source Code may also be made available under the following Secondary
 * Licenses when the conditions for such availability set forth in the Eclipse
 * Public License v. 2.0 are satisfied: GPL-3.0 which is available at
 * https://www.gnu.org/licenses/gpl-3.0-standalone.html.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
grammar InternalElkGraphJson;

options {
	superClass=AbstractInternalContentAssistParser;
}

@lexer::header {
package org.eclipse.elk.graph.json.text.ide.contentassist.antlr.internal;

// Hack: Use our own Lexer superclass by means of import. 
// Currently there is no other way to specify the superclass for the lexer.
import org.eclipse.xtext.ide.editor.contentassist.antlr.internal.Lexer;
}

@parser::header {
package org.eclipse.elk.graph.json.text.ide.contentassist.antlr.internal;

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
import org.eclipse.elk.graph.json.text.services.ElkGraphJsonGrammarAccess;

}
@parser::members {
	private ElkGraphJsonGrammarAccess grammarAccess;

	public void setGrammarAccess(ElkGraphJsonGrammarAccess grammarAccess) {
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


// Rule NodeElement
ruleNodeElement 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getNodeElementAccess().getAlternatives()); }
		(rule__NodeElement__Alternatives)
		{ after(grammarAccess.getNodeElementAccess().getAlternatives()); }
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


// Rule PortElement
rulePortElement 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getPortElementAccess().getAlternatives()); }
		(rule__PortElement__Alternatives)
		{ after(grammarAccess.getPortElementAccess().getAlternatives()); }
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


// Rule LabelElement
ruleLabelElement 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getLabelElementAccess().getAlternatives()); }
		(rule__LabelElement__Alternatives)
		{ after(grammarAccess.getLabelElementAccess().getAlternatives()); }
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


// Rule EdgeElement
ruleEdgeElement 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getEdgeElementAccess().getAlternatives()); }
		(rule__EdgeElement__Alternatives)
		{ after(grammarAccess.getEdgeElementAccess().getAlternatives()); }
	)
;
finally {
	restoreStackSize(stackSize);
}


// Rule ElkEdgeSources
ruleElkEdgeSources 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getElkEdgeSourcesAccess().getGroup()); }
		(rule__ElkEdgeSources__Group__0)
		{ after(grammarAccess.getElkEdgeSourcesAccess().getGroup()); }
	)
;
finally {
	restoreStackSize(stackSize);
}


// Rule ElkEdgeTargets
ruleElkEdgeTargets 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getElkEdgeTargetsAccess().getGroup()); }
		(rule__ElkEdgeTargets__Group__0)
		{ after(grammarAccess.getElkEdgeTargetsAccess().getGroup()); }
	)
;
finally {
	restoreStackSize(stackSize);
}


// Rule ElkId
ruleElkId 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getElkIdAccess().getGroup()); }
		(rule__ElkId__Group__0)
		{ after(grammarAccess.getElkIdAccess().getGroup()); }
	)
;
finally {
	restoreStackSize(stackSize);
}


// Rule ElkNodeChildren
ruleElkNodeChildren 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getElkNodeChildrenAccess().getGroup()); }
		(rule__ElkNodeChildren__Group__0)
		{ after(grammarAccess.getElkNodeChildrenAccess().getGroup()); }
	)
;
finally {
	restoreStackSize(stackSize);
}


// Rule ElkNodePorts
ruleElkNodePorts 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getElkNodePortsAccess().getGroup()); }
		(rule__ElkNodePorts__Group__0)
		{ after(grammarAccess.getElkNodePortsAccess().getGroup()); }
	)
;
finally {
	restoreStackSize(stackSize);
}


// Rule ElkNodeEdges
ruleElkNodeEdges 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getElkNodeEdgesAccess().getGroup()); }
		(rule__ElkNodeEdges__Group__0)
		{ after(grammarAccess.getElkNodeEdgesAccess().getGroup()); }
	)
;
finally {
	restoreStackSize(stackSize);
}


// Rule ElkGraphElementLabels
ruleElkGraphElementLabels 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getElkGraphElementLabelsAccess().getGroup()); }
		(rule__ElkGraphElementLabels__Group__0)
		{ after(grammarAccess.getElkGraphElementLabelsAccess().getGroup()); }
	)
;
finally {
	restoreStackSize(stackSize);
}


// Rule ElkGraphElementProperties
ruleElkGraphElementProperties 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getElkGraphElementPropertiesAccess().getGroup()); }
		(rule__ElkGraphElementProperties__Group__0)
		{ after(grammarAccess.getElkGraphElementPropertiesAccess().getGroup()); }
	)
;
finally {
	restoreStackSize(stackSize);
}


// Rule ShapeElement
ruleShapeElement 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getShapeElementAccess().getAlternatives()); }
		(rule__ShapeElement__Alternatives)
		{ after(grammarAccess.getShapeElementAccess().getAlternatives()); }
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
		{ before(grammarAccess.getPropertyKeyAccess().getAlternatives()); }
		(rule__PropertyKey__Alternatives)
		{ after(grammarAccess.getPropertyKeyAccess().getAlternatives()); }
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

// Entry rule entryRuleJsonObject
entryRuleJsonObject
:
{ before(grammarAccess.getJsonObjectRule()); }
	 ruleJsonObject
{ after(grammarAccess.getJsonObjectRule()); } 
	 EOF 
;

// Rule JsonObject
ruleJsonObject 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getJsonObjectAccess().getGroup()); }
		(rule__JsonObject__Group__0)
		{ after(grammarAccess.getJsonObjectAccess().getGroup()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

// Entry rule entryRuleJsonArray
entryRuleJsonArray
:
{ before(grammarAccess.getJsonArrayRule()); }
	 ruleJsonArray
{ after(grammarAccess.getJsonArrayRule()); } 
	 EOF 
;

// Rule JsonArray
ruleJsonArray 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getJsonArrayAccess().getGroup()); }
		(rule__JsonArray__Group__0)
		{ after(grammarAccess.getJsonArrayAccess().getGroup()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

// Entry rule entryRuleJsonMember
entryRuleJsonMember
:
{ before(grammarAccess.getJsonMemberRule()); }
	 ruleJsonMember
{ after(grammarAccess.getJsonMemberRule()); } 
	 EOF 
;

// Rule JsonMember
ruleJsonMember 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getJsonMemberAccess().getGroup()); }
		(rule__JsonMember__Group__0)
		{ after(grammarAccess.getJsonMemberAccess().getGroup()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

// Entry rule entryRuleJsonValue
entryRuleJsonValue
:
{ before(grammarAccess.getJsonValueRule()); }
	 ruleJsonValue
{ after(grammarAccess.getJsonValueRule()); } 
	 EOF 
;

// Rule JsonValue
ruleJsonValue 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getJsonValueAccess().getAlternatives()); }
		(rule__JsonValue__Alternatives)
		{ after(grammarAccess.getJsonValueAccess().getAlternatives()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

// Entry rule entryRuleKeyChildren
entryRuleKeyChildren
:
{ before(grammarAccess.getKeyChildrenRule()); }
	 ruleKeyChildren
{ after(grammarAccess.getKeyChildrenRule()); } 
	 EOF 
;

// Rule KeyChildren
ruleKeyChildren 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getKeyChildrenAccess().getAlternatives()); }
		(rule__KeyChildren__Alternatives)
		{ after(grammarAccess.getKeyChildrenAccess().getAlternatives()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

// Entry rule entryRuleKeyPorts
entryRuleKeyPorts
:
{ before(grammarAccess.getKeyPortsRule()); }
	 ruleKeyPorts
{ after(grammarAccess.getKeyPortsRule()); } 
	 EOF 
;

// Rule KeyPorts
ruleKeyPorts 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getKeyPortsAccess().getAlternatives()); }
		(rule__KeyPorts__Alternatives)
		{ after(grammarAccess.getKeyPortsAccess().getAlternatives()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

// Entry rule entryRuleKeyLabels
entryRuleKeyLabels
:
{ before(grammarAccess.getKeyLabelsRule()); }
	 ruleKeyLabels
{ after(grammarAccess.getKeyLabelsRule()); } 
	 EOF 
;

// Rule KeyLabels
ruleKeyLabels 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getKeyLabelsAccess().getAlternatives()); }
		(rule__KeyLabels__Alternatives)
		{ after(grammarAccess.getKeyLabelsAccess().getAlternatives()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

// Entry rule entryRuleKeyEdges
entryRuleKeyEdges
:
{ before(grammarAccess.getKeyEdgesRule()); }
	 ruleKeyEdges
{ after(grammarAccess.getKeyEdgesRule()); } 
	 EOF 
;

// Rule KeyEdges
ruleKeyEdges 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getKeyEdgesAccess().getAlternatives()); }
		(rule__KeyEdges__Alternatives)
		{ after(grammarAccess.getKeyEdgesAccess().getAlternatives()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

// Entry rule entryRuleKeyLayoutOptions
entryRuleKeyLayoutOptions
:
{ before(grammarAccess.getKeyLayoutOptionsRule()); }
	 ruleKeyLayoutOptions
{ after(grammarAccess.getKeyLayoutOptionsRule()); } 
	 EOF 
;

// Rule KeyLayoutOptions
ruleKeyLayoutOptions 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getKeyLayoutOptionsAccess().getAlternatives()); }
		(rule__KeyLayoutOptions__Alternatives)
		{ after(grammarAccess.getKeyLayoutOptionsAccess().getAlternatives()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

// Entry rule entryRuleKeyId
entryRuleKeyId
:
{ before(grammarAccess.getKeyIdRule()); }
	 ruleKeyId
{ after(grammarAccess.getKeyIdRule()); } 
	 EOF 
;

// Rule KeyId
ruleKeyId 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getKeyIdAccess().getAlternatives()); }
		(rule__KeyId__Alternatives)
		{ after(grammarAccess.getKeyIdAccess().getAlternatives()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

// Entry rule entryRuleKeyX
entryRuleKeyX
:
{ before(grammarAccess.getKeyXRule()); }
	 ruleKeyX
{ after(grammarAccess.getKeyXRule()); } 
	 EOF 
;

// Rule KeyX
ruleKeyX 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getKeyXAccess().getAlternatives()); }
		(rule__KeyX__Alternatives)
		{ after(grammarAccess.getKeyXAccess().getAlternatives()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

// Entry rule entryRuleKeyY
entryRuleKeyY
:
{ before(grammarAccess.getKeyYRule()); }
	 ruleKeyY
{ after(grammarAccess.getKeyYRule()); } 
	 EOF 
;

// Rule KeyY
ruleKeyY 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getKeyYAccess().getAlternatives()); }
		(rule__KeyY__Alternatives)
		{ after(grammarAccess.getKeyYAccess().getAlternatives()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

// Entry rule entryRuleKeyWidth
entryRuleKeyWidth
:
{ before(grammarAccess.getKeyWidthRule()); }
	 ruleKeyWidth
{ after(grammarAccess.getKeyWidthRule()); } 
	 EOF 
;

// Rule KeyWidth
ruleKeyWidth 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getKeyWidthAccess().getAlternatives()); }
		(rule__KeyWidth__Alternatives)
		{ after(grammarAccess.getKeyWidthAccess().getAlternatives()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

// Entry rule entryRuleKeyHeight
entryRuleKeyHeight
:
{ before(grammarAccess.getKeyHeightRule()); }
	 ruleKeyHeight
{ after(grammarAccess.getKeyHeightRule()); } 
	 EOF 
;

// Rule KeyHeight
ruleKeyHeight 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getKeyHeightAccess().getAlternatives()); }
		(rule__KeyHeight__Alternatives)
		{ after(grammarAccess.getKeyHeightAccess().getAlternatives()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

// Entry rule entryRuleKeySources
entryRuleKeySources
:
{ before(grammarAccess.getKeySourcesRule()); }
	 ruleKeySources
{ after(grammarAccess.getKeySourcesRule()); } 
	 EOF 
;

// Rule KeySources
ruleKeySources 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getKeySourcesAccess().getAlternatives()); }
		(rule__KeySources__Alternatives)
		{ after(grammarAccess.getKeySourcesAccess().getAlternatives()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

// Entry rule entryRuleKeyTargets
entryRuleKeyTargets
:
{ before(grammarAccess.getKeyTargetsRule()); }
	 ruleKeyTargets
{ after(grammarAccess.getKeyTargetsRule()); } 
	 EOF 
;

// Rule KeyTargets
ruleKeyTargets 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getKeyTargetsAccess().getAlternatives()); }
		(rule__KeyTargets__Alternatives)
		{ after(grammarAccess.getKeyTargetsAccess().getAlternatives()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

// Entry rule entryRuleKeyText
entryRuleKeyText
:
{ before(grammarAccess.getKeyTextRule()); }
	 ruleKeyText
{ after(grammarAccess.getKeyTextRule()); } 
	 EOF 
;

// Rule KeyText
ruleKeyText 
	@init {
		int stackSize = keepStackSize();
	}
	:
	(
		{ before(grammarAccess.getKeyTextAccess().getAlternatives()); }
		(rule__KeyText__Alternatives)
		{ after(grammarAccess.getKeyTextAccess().getAlternatives()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__NodeElement__Alternatives
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getNodeElementAccess().getElkIdParserRuleCall_0()); }
		ruleElkId
		{ after(grammarAccess.getNodeElementAccess().getElkIdParserRuleCall_0()); }
	)
	|
	(
		{ before(grammarAccess.getNodeElementAccess().getShapeElementParserRuleCall_1()); }
		ruleShapeElement
		{ after(grammarAccess.getNodeElementAccess().getShapeElementParserRuleCall_1()); }
	)
	|
	(
		{ before(grammarAccess.getNodeElementAccess().getGroup_2()); }
		(rule__NodeElement__Group_2__0)
		{ after(grammarAccess.getNodeElementAccess().getGroup_2()); }
	)
	|
	(
		{ before(grammarAccess.getNodeElementAccess().getGroup_3()); }
		(rule__NodeElement__Group_3__0)
		{ after(grammarAccess.getNodeElementAccess().getGroup_3()); }
	)
	|
	(
		{ before(grammarAccess.getNodeElementAccess().getGroup_4()); }
		(rule__NodeElement__Group_4__0)
		{ after(grammarAccess.getNodeElementAccess().getGroup_4()); }
	)
	|
	(
		{ before(grammarAccess.getNodeElementAccess().getGroup_5()); }
		(rule__NodeElement__Group_5__0)
		{ after(grammarAccess.getNodeElementAccess().getGroup_5()); }
	)
	|
	(
		{ before(grammarAccess.getNodeElementAccess().getGroup_6()); }
		(rule__NodeElement__Group_6__0)
		{ after(grammarAccess.getNodeElementAccess().getGroup_6()); }
	)
	|
	(
		{ before(grammarAccess.getNodeElementAccess().getJsonMemberParserRuleCall_7()); }
		ruleJsonMember
		{ after(grammarAccess.getNodeElementAccess().getJsonMemberParserRuleCall_7()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__PortElement__Alternatives
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getPortElementAccess().getElkIdParserRuleCall_0()); }
		ruleElkId
		{ after(grammarAccess.getPortElementAccess().getElkIdParserRuleCall_0()); }
	)
	|
	(
		{ before(grammarAccess.getPortElementAccess().getShapeElementParserRuleCall_1()); }
		ruleShapeElement
		{ after(grammarAccess.getPortElementAccess().getShapeElementParserRuleCall_1()); }
	)
	|
	(
		{ before(grammarAccess.getPortElementAccess().getGroup_2()); }
		(rule__PortElement__Group_2__0)
		{ after(grammarAccess.getPortElementAccess().getGroup_2()); }
	)
	|
	(
		{ before(grammarAccess.getPortElementAccess().getGroup_3()); }
		(rule__PortElement__Group_3__0)
		{ after(grammarAccess.getPortElementAccess().getGroup_3()); }
	)
	|
	(
		{ before(grammarAccess.getPortElementAccess().getJsonMemberParserRuleCall_4()); }
		ruleJsonMember
		{ after(grammarAccess.getPortElementAccess().getJsonMemberParserRuleCall_4()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__LabelElement__Alternatives
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getLabelElementAccess().getElkIdParserRuleCall_0()); }
		ruleElkId
		{ after(grammarAccess.getLabelElementAccess().getElkIdParserRuleCall_0()); }
	)
	|
	(
		{ before(grammarAccess.getLabelElementAccess().getShapeElementParserRuleCall_1()); }
		ruleShapeElement
		{ after(grammarAccess.getLabelElementAccess().getShapeElementParserRuleCall_1()); }
	)
	|
	(
		{ before(grammarAccess.getLabelElementAccess().getGroup_2()); }
		(rule__LabelElement__Group_2__0)
		{ after(grammarAccess.getLabelElementAccess().getGroup_2()); }
	)
	|
	(
		{ before(grammarAccess.getLabelElementAccess().getGroup_3()); }
		(rule__LabelElement__Group_3__0)
		{ after(grammarAccess.getLabelElementAccess().getGroup_3()); }
	)
	|
	(
		{ before(grammarAccess.getLabelElementAccess().getGroup_4()); }
		(rule__LabelElement__Group_4__0)
		{ after(grammarAccess.getLabelElementAccess().getGroup_4()); }
	)
	|
	(
		{ before(grammarAccess.getLabelElementAccess().getJsonMemberParserRuleCall_5()); }
		ruleJsonMember
		{ after(grammarAccess.getLabelElementAccess().getJsonMemberParserRuleCall_5()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__EdgeElement__Alternatives
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getEdgeElementAccess().getElkIdParserRuleCall_0()); }
		ruleElkId
		{ after(grammarAccess.getEdgeElementAccess().getElkIdParserRuleCall_0()); }
	)
	|
	(
		{ before(grammarAccess.getEdgeElementAccess().getGroup_1()); }
		(rule__EdgeElement__Group_1__0)
		{ after(grammarAccess.getEdgeElementAccess().getGroup_1()); }
	)
	|
	(
		{ before(grammarAccess.getEdgeElementAccess().getGroup_2()); }
		(rule__EdgeElement__Group_2__0)
		{ after(grammarAccess.getEdgeElementAccess().getGroup_2()); }
	)
	|
	(
		{ before(grammarAccess.getEdgeElementAccess().getGroup_3()); }
		(rule__EdgeElement__Group_3__0)
		{ after(grammarAccess.getEdgeElementAccess().getGroup_3()); }
	)
	|
	(
		{ before(grammarAccess.getEdgeElementAccess().getGroup_4()); }
		(rule__EdgeElement__Group_4__0)
		{ after(grammarAccess.getEdgeElementAccess().getGroup_4()); }
	)
	|
	(
		{ before(grammarAccess.getEdgeElementAccess().getJsonMemberParserRuleCall_5()); }
		ruleJsonMember
		{ after(grammarAccess.getEdgeElementAccess().getJsonMemberParserRuleCall_5()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeElement__Alternatives
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getShapeElementAccess().getGroup_0()); }
		(rule__ShapeElement__Group_0__0)
		{ after(grammarAccess.getShapeElementAccess().getGroup_0()); }
	)
	|
	(
		{ before(grammarAccess.getShapeElementAccess().getGroup_1()); }
		(rule__ShapeElement__Group_1__0)
		{ after(grammarAccess.getShapeElementAccess().getGroup_1()); }
	)
	|
	(
		{ before(grammarAccess.getShapeElementAccess().getGroup_2()); }
		(rule__ShapeElement__Group_2__0)
		{ after(grammarAccess.getShapeElementAccess().getGroup_2()); }
	)
	|
	(
		{ before(grammarAccess.getShapeElementAccess().getGroup_3()); }
		(rule__ShapeElement__Group_3__0)
		{ after(grammarAccess.getShapeElementAccess().getGroup_3()); }
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
;
finally {
	restoreStackSize(stackSize);
}

rule__PropertyKey__Alternatives
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getPropertyKeyAccess().getSTRINGTerminalRuleCall_0()); }
		RULE_STRING
		{ after(grammarAccess.getPropertyKeyAccess().getSTRINGTerminalRuleCall_0()); }
	)
	|
	(
		{ before(grammarAccess.getPropertyKeyAccess().getIDTerminalRuleCall_1()); }
		RULE_ID
		{ after(grammarAccess.getPropertyKeyAccess().getIDTerminalRuleCall_1()); }
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

rule__JsonMember__Alternatives_0
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getJsonMemberAccess().getSTRINGTerminalRuleCall_0_0()); }
		RULE_STRING
		{ after(grammarAccess.getJsonMemberAccess().getSTRINGTerminalRuleCall_0_0()); }
	)
	|
	(
		{ before(grammarAccess.getJsonMemberAccess().getIDTerminalRuleCall_0_1()); }
		RULE_ID
		{ after(grammarAccess.getJsonMemberAccess().getIDTerminalRuleCall_0_1()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__JsonValue__Alternatives
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getJsonValueAccess().getSTRINGTerminalRuleCall_0()); }
		RULE_STRING
		{ after(grammarAccess.getJsonValueAccess().getSTRINGTerminalRuleCall_0()); }
	)
	|
	(
		{ before(grammarAccess.getJsonValueAccess().getNumberParserRuleCall_1()); }
		ruleNumber
		{ after(grammarAccess.getJsonValueAccess().getNumberParserRuleCall_1()); }
	)
	|
	(
		{ before(grammarAccess.getJsonValueAccess().getJsonObjectParserRuleCall_2()); }
		ruleJsonObject
		{ after(grammarAccess.getJsonValueAccess().getJsonObjectParserRuleCall_2()); }
	)
	|
	(
		{ before(grammarAccess.getJsonValueAccess().getJsonArrayParserRuleCall_3()); }
		ruleJsonArray
		{ after(grammarAccess.getJsonValueAccess().getJsonArrayParserRuleCall_3()); }
	)
	|
	(
		{ before(grammarAccess.getJsonValueAccess().getBooleanValueParserRuleCall_4()); }
		ruleBooleanValue
		{ after(grammarAccess.getJsonValueAccess().getBooleanValueParserRuleCall_4()); }
	)
	|
	(
		{ before(grammarAccess.getJsonValueAccess().getNullKeyword_5()); }
		'null'
		{ after(grammarAccess.getJsonValueAccess().getNullKeyword_5()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__KeyChildren__Alternatives
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getKeyChildrenAccess().getChildrenKeyword_0()); }
		'"children"'
		{ after(grammarAccess.getKeyChildrenAccess().getChildrenKeyword_0()); }
	)
	|
	(
		{ before(grammarAccess.getKeyChildrenAccess().getChildrenKeyword_1()); }
		'\'children\''
		{ after(grammarAccess.getKeyChildrenAccess().getChildrenKeyword_1()); }
	)
	|
	(
		{ before(grammarAccess.getKeyChildrenAccess().getChildrenKeyword_2()); }
		'children'
		{ after(grammarAccess.getKeyChildrenAccess().getChildrenKeyword_2()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__KeyPorts__Alternatives
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getKeyPortsAccess().getPortsKeyword_0()); }
		'"ports"'
		{ after(grammarAccess.getKeyPortsAccess().getPortsKeyword_0()); }
	)
	|
	(
		{ before(grammarAccess.getKeyPortsAccess().getPortsKeyword_1()); }
		'\'ports\''
		{ after(grammarAccess.getKeyPortsAccess().getPortsKeyword_1()); }
	)
	|
	(
		{ before(grammarAccess.getKeyPortsAccess().getPortsKeyword_2()); }
		'ports'
		{ after(grammarAccess.getKeyPortsAccess().getPortsKeyword_2()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__KeyLabels__Alternatives
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getKeyLabelsAccess().getLabelsKeyword_0()); }
		'"labels"'
		{ after(grammarAccess.getKeyLabelsAccess().getLabelsKeyword_0()); }
	)
	|
	(
		{ before(grammarAccess.getKeyLabelsAccess().getLabelsKeyword_1()); }
		'\'labels\''
		{ after(grammarAccess.getKeyLabelsAccess().getLabelsKeyword_1()); }
	)
	|
	(
		{ before(grammarAccess.getKeyLabelsAccess().getLabelsKeyword_2()); }
		'labels'
		{ after(grammarAccess.getKeyLabelsAccess().getLabelsKeyword_2()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__KeyEdges__Alternatives
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getKeyEdgesAccess().getEdgesKeyword_0()); }
		'"edges"'
		{ after(grammarAccess.getKeyEdgesAccess().getEdgesKeyword_0()); }
	)
	|
	(
		{ before(grammarAccess.getKeyEdgesAccess().getEdgesKeyword_1()); }
		'\'edges\''
		{ after(grammarAccess.getKeyEdgesAccess().getEdgesKeyword_1()); }
	)
	|
	(
		{ before(grammarAccess.getKeyEdgesAccess().getEdgesKeyword_2()); }
		'edges'
		{ after(grammarAccess.getKeyEdgesAccess().getEdgesKeyword_2()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__KeyLayoutOptions__Alternatives
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getKeyLayoutOptionsAccess().getLayoutOptionsKeyword_0()); }
		'"layoutOptions"'
		{ after(grammarAccess.getKeyLayoutOptionsAccess().getLayoutOptionsKeyword_0()); }
	)
	|
	(
		{ before(grammarAccess.getKeyLayoutOptionsAccess().getLayoutOptionsKeyword_1()); }
		'\'layoutOptions\''
		{ after(grammarAccess.getKeyLayoutOptionsAccess().getLayoutOptionsKeyword_1()); }
	)
	|
	(
		{ before(grammarAccess.getKeyLayoutOptionsAccess().getLayoutOptionsKeyword_2()); }
		'layoutOptions'
		{ after(grammarAccess.getKeyLayoutOptionsAccess().getLayoutOptionsKeyword_2()); }
	)
	|
	(
		{ before(grammarAccess.getKeyLayoutOptionsAccess().getPropertiesKeyword_3()); }
		'"properties"'
		{ after(grammarAccess.getKeyLayoutOptionsAccess().getPropertiesKeyword_3()); }
	)
	|
	(
		{ before(grammarAccess.getKeyLayoutOptionsAccess().getPropertiesKeyword_4()); }
		'\'properties\''
		{ after(grammarAccess.getKeyLayoutOptionsAccess().getPropertiesKeyword_4()); }
	)
	|
	(
		{ before(grammarAccess.getKeyLayoutOptionsAccess().getPropertiesKeyword_5()); }
		'properties'
		{ after(grammarAccess.getKeyLayoutOptionsAccess().getPropertiesKeyword_5()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__KeyId__Alternatives
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getKeyIdAccess().getIdKeyword_0()); }
		'"id"'
		{ after(grammarAccess.getKeyIdAccess().getIdKeyword_0()); }
	)
	|
	(
		{ before(grammarAccess.getKeyIdAccess().getIdKeyword_1()); }
		'\'id\''
		{ after(grammarAccess.getKeyIdAccess().getIdKeyword_1()); }
	)
	|
	(
		{ before(grammarAccess.getKeyIdAccess().getIdKeyword_2()); }
		'id'
		{ after(grammarAccess.getKeyIdAccess().getIdKeyword_2()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__KeyX__Alternatives
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getKeyXAccess().getXKeyword_0()); }
		'"x"'
		{ after(grammarAccess.getKeyXAccess().getXKeyword_0()); }
	)
	|
	(
		{ before(grammarAccess.getKeyXAccess().getXKeyword_1()); }
		'\'x\''
		{ after(grammarAccess.getKeyXAccess().getXKeyword_1()); }
	)
	|
	(
		{ before(grammarAccess.getKeyXAccess().getXKeyword_2()); }
		'x'
		{ after(grammarAccess.getKeyXAccess().getXKeyword_2()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__KeyY__Alternatives
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getKeyYAccess().getYKeyword_0()); }
		'"y"'
		{ after(grammarAccess.getKeyYAccess().getYKeyword_0()); }
	)
	|
	(
		{ before(grammarAccess.getKeyYAccess().getYKeyword_1()); }
		'\'y\''
		{ after(grammarAccess.getKeyYAccess().getYKeyword_1()); }
	)
	|
	(
		{ before(grammarAccess.getKeyYAccess().getYKeyword_2()); }
		'y'
		{ after(grammarAccess.getKeyYAccess().getYKeyword_2()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__KeyWidth__Alternatives
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getKeyWidthAccess().getWidthKeyword_0()); }
		'"width"'
		{ after(grammarAccess.getKeyWidthAccess().getWidthKeyword_0()); }
	)
	|
	(
		{ before(grammarAccess.getKeyWidthAccess().getWidthKeyword_1()); }
		'\'width\''
		{ after(grammarAccess.getKeyWidthAccess().getWidthKeyword_1()); }
	)
	|
	(
		{ before(grammarAccess.getKeyWidthAccess().getWidthKeyword_2()); }
		'width'
		{ after(grammarAccess.getKeyWidthAccess().getWidthKeyword_2()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__KeyHeight__Alternatives
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getKeyHeightAccess().getHeightKeyword_0()); }
		'"height"'
		{ after(grammarAccess.getKeyHeightAccess().getHeightKeyword_0()); }
	)
	|
	(
		{ before(grammarAccess.getKeyHeightAccess().getHeightKeyword_1()); }
		'\'height\''
		{ after(grammarAccess.getKeyHeightAccess().getHeightKeyword_1()); }
	)
	|
	(
		{ before(grammarAccess.getKeyHeightAccess().getHeightKeyword_2()); }
		'height'
		{ after(grammarAccess.getKeyHeightAccess().getHeightKeyword_2()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__KeySources__Alternatives
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getKeySourcesAccess().getSourcesKeyword_0()); }
		'"sources"'
		{ after(grammarAccess.getKeySourcesAccess().getSourcesKeyword_0()); }
	)
	|
	(
		{ before(grammarAccess.getKeySourcesAccess().getSourcesKeyword_1()); }
		'\'sources\''
		{ after(grammarAccess.getKeySourcesAccess().getSourcesKeyword_1()); }
	)
	|
	(
		{ before(grammarAccess.getKeySourcesAccess().getSourcesKeyword_2()); }
		'sources'
		{ after(grammarAccess.getKeySourcesAccess().getSourcesKeyword_2()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__KeyTargets__Alternatives
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getKeyTargetsAccess().getTargetsKeyword_0()); }
		'"targets"'
		{ after(grammarAccess.getKeyTargetsAccess().getTargetsKeyword_0()); }
	)
	|
	(
		{ before(grammarAccess.getKeyTargetsAccess().getTargetsKeyword_1()); }
		'\'targets\''
		{ after(grammarAccess.getKeyTargetsAccess().getTargetsKeyword_1()); }
	)
	|
	(
		{ before(grammarAccess.getKeyTargetsAccess().getTargetsKeyword_2()); }
		'targets'
		{ after(grammarAccess.getKeyTargetsAccess().getTargetsKeyword_2()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__KeyText__Alternatives
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getKeyTextAccess().getTextKeyword_0()); }
		'"text"'
		{ after(grammarAccess.getKeyTextAccess().getTextKeyword_0()); }
	)
	|
	(
		{ before(grammarAccess.getKeyTextAccess().getTextKeyword_1()); }
		'\'text\''
		{ after(grammarAccess.getKeyTextAccess().getTextKeyword_1()); }
	)
	|
	(
		{ before(grammarAccess.getKeyTextAccess().getTextKeyword_2()); }
		'text'
		{ after(grammarAccess.getKeyTextAccess().getTextKeyword_2()); }
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
	{ before(grammarAccess.getElkNodeAccess().getElkNodeAction_0()); }
	()
	{ after(grammarAccess.getElkNodeAccess().getElkNodeAction_0()); }
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
	{ before(grammarAccess.getElkNodeAccess().getLeftCurlyBracketKeyword_1()); }
	'{'
	{ after(grammarAccess.getElkNodeAccess().getLeftCurlyBracketKeyword_1()); }
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
	rule__ElkNode__Group__3
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

rule__ElkNode__Group__3
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkNode__Group__3__Impl
	rule__ElkNode__Group__4
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNode__Group__3__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkNodeAccess().getCommaKeyword_3()); }
	(',')?
	{ after(grammarAccess.getElkNodeAccess().getCommaKeyword_3()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNode__Group__4
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkNode__Group__4__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNode__Group__4__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkNodeAccess().getRightCurlyBracketKeyword_4()); }
	'}'
	{ after(grammarAccess.getElkNodeAccess().getRightCurlyBracketKeyword_4()); }
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
	{ before(grammarAccess.getElkNodeAccess().getNodeElementParserRuleCall_2_0()); }
	ruleNodeElement
	{ after(grammarAccess.getElkNodeAccess().getNodeElementParserRuleCall_2_0()); }
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
	{ before(grammarAccess.getElkNodeAccess().getGroup_2_1()); }
	(rule__ElkNode__Group_2_1__0)*
	{ after(grammarAccess.getElkNodeAccess().getGroup_2_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkNode__Group_2_1__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkNode__Group_2_1__0__Impl
	rule__ElkNode__Group_2_1__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNode__Group_2_1__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkNodeAccess().getCommaKeyword_2_1_0()); }
	','
	{ after(grammarAccess.getElkNodeAccess().getCommaKeyword_2_1_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNode__Group_2_1__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkNode__Group_2_1__1__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNode__Group_2_1__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkNodeAccess().getNodeElementParserRuleCall_2_1_1()); }
	ruleNodeElement
	{ after(grammarAccess.getElkNodeAccess().getNodeElementParserRuleCall_2_1_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__NodeElement__Group_2__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__NodeElement__Group_2__0__Impl
	rule__NodeElement__Group_2__1
;
finally {
	restoreStackSize(stackSize);
}

rule__NodeElement__Group_2__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getNodeElementAccess().getKeyChildrenParserRuleCall_2_0()); }
	ruleKeyChildren
	{ after(grammarAccess.getNodeElementAccess().getKeyChildrenParserRuleCall_2_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__NodeElement__Group_2__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__NodeElement__Group_2__1__Impl
	rule__NodeElement__Group_2__2
;
finally {
	restoreStackSize(stackSize);
}

rule__NodeElement__Group_2__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getNodeElementAccess().getColonKeyword_2_1()); }
	':'
	{ after(grammarAccess.getNodeElementAccess().getColonKeyword_2_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__NodeElement__Group_2__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__NodeElement__Group_2__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__NodeElement__Group_2__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getNodeElementAccess().getElkNodeChildrenParserRuleCall_2_2()); }
	ruleElkNodeChildren
	{ after(grammarAccess.getNodeElementAccess().getElkNodeChildrenParserRuleCall_2_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__NodeElement__Group_3__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__NodeElement__Group_3__0__Impl
	rule__NodeElement__Group_3__1
;
finally {
	restoreStackSize(stackSize);
}

rule__NodeElement__Group_3__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getNodeElementAccess().getKeyPortsParserRuleCall_3_0()); }
	ruleKeyPorts
	{ after(grammarAccess.getNodeElementAccess().getKeyPortsParserRuleCall_3_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__NodeElement__Group_3__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__NodeElement__Group_3__1__Impl
	rule__NodeElement__Group_3__2
;
finally {
	restoreStackSize(stackSize);
}

rule__NodeElement__Group_3__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getNodeElementAccess().getColonKeyword_3_1()); }
	':'
	{ after(grammarAccess.getNodeElementAccess().getColonKeyword_3_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__NodeElement__Group_3__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__NodeElement__Group_3__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__NodeElement__Group_3__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getNodeElementAccess().getElkNodePortsParserRuleCall_3_2()); }
	ruleElkNodePorts
	{ after(grammarAccess.getNodeElementAccess().getElkNodePortsParserRuleCall_3_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__NodeElement__Group_4__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__NodeElement__Group_4__0__Impl
	rule__NodeElement__Group_4__1
;
finally {
	restoreStackSize(stackSize);
}

rule__NodeElement__Group_4__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getNodeElementAccess().getKeyLabelsParserRuleCall_4_0()); }
	ruleKeyLabels
	{ after(grammarAccess.getNodeElementAccess().getKeyLabelsParserRuleCall_4_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__NodeElement__Group_4__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__NodeElement__Group_4__1__Impl
	rule__NodeElement__Group_4__2
;
finally {
	restoreStackSize(stackSize);
}

rule__NodeElement__Group_4__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getNodeElementAccess().getColonKeyword_4_1()); }
	':'
	{ after(grammarAccess.getNodeElementAccess().getColonKeyword_4_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__NodeElement__Group_4__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__NodeElement__Group_4__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__NodeElement__Group_4__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getNodeElementAccess().getElkGraphElementLabelsParserRuleCall_4_2()); }
	ruleElkGraphElementLabels
	{ after(grammarAccess.getNodeElementAccess().getElkGraphElementLabelsParserRuleCall_4_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__NodeElement__Group_5__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__NodeElement__Group_5__0__Impl
	rule__NodeElement__Group_5__1
;
finally {
	restoreStackSize(stackSize);
}

rule__NodeElement__Group_5__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getNodeElementAccess().getKeyEdgesParserRuleCall_5_0()); }
	ruleKeyEdges
	{ after(grammarAccess.getNodeElementAccess().getKeyEdgesParserRuleCall_5_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__NodeElement__Group_5__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__NodeElement__Group_5__1__Impl
	rule__NodeElement__Group_5__2
;
finally {
	restoreStackSize(stackSize);
}

rule__NodeElement__Group_5__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getNodeElementAccess().getColonKeyword_5_1()); }
	':'
	{ after(grammarAccess.getNodeElementAccess().getColonKeyword_5_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__NodeElement__Group_5__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__NodeElement__Group_5__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__NodeElement__Group_5__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getNodeElementAccess().getElkNodeEdgesParserRuleCall_5_2()); }
	ruleElkNodeEdges
	{ after(grammarAccess.getNodeElementAccess().getElkNodeEdgesParserRuleCall_5_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__NodeElement__Group_6__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__NodeElement__Group_6__0__Impl
	rule__NodeElement__Group_6__1
;
finally {
	restoreStackSize(stackSize);
}

rule__NodeElement__Group_6__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getNodeElementAccess().getKeyLayoutOptionsParserRuleCall_6_0()); }
	ruleKeyLayoutOptions
	{ after(grammarAccess.getNodeElementAccess().getKeyLayoutOptionsParserRuleCall_6_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__NodeElement__Group_6__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__NodeElement__Group_6__1__Impl
	rule__NodeElement__Group_6__2
;
finally {
	restoreStackSize(stackSize);
}

rule__NodeElement__Group_6__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getNodeElementAccess().getColonKeyword_6_1()); }
	':'
	{ after(grammarAccess.getNodeElementAccess().getColonKeyword_6_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__NodeElement__Group_6__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__NodeElement__Group_6__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__NodeElement__Group_6__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getNodeElementAccess().getElkGraphElementPropertiesParserRuleCall_6_2()); }
	ruleElkGraphElementProperties
	{ after(grammarAccess.getNodeElementAccess().getElkGraphElementPropertiesParserRuleCall_6_2()); }
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
	{ before(grammarAccess.getElkPortAccess().getLeftCurlyBracketKeyword_0()); }
	'{'
	{ after(grammarAccess.getElkPortAccess().getLeftCurlyBracketKeyword_0()); }
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
	{ before(grammarAccess.getElkPortAccess().getPortElementParserRuleCall_1()); }
	rulePortElement
	{ after(grammarAccess.getElkPortAccess().getPortElementParserRuleCall_1()); }
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
	rule__ElkPort__Group__3
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
	(rule__ElkPort__Group_2__0)*
	{ after(grammarAccess.getElkPortAccess().getGroup_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkPort__Group__3
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkPort__Group__3__Impl
	rule__ElkPort__Group__4
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkPort__Group__3__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkPortAccess().getCommaKeyword_3()); }
	(',')?
	{ after(grammarAccess.getElkPortAccess().getCommaKeyword_3()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkPort__Group__4
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkPort__Group__4__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkPort__Group__4__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkPortAccess().getRightCurlyBracketKeyword_4()); }
	'}'
	{ after(grammarAccess.getElkPortAccess().getRightCurlyBracketKeyword_4()); }
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
	{ before(grammarAccess.getElkPortAccess().getCommaKeyword_2_0()); }
	','
	{ after(grammarAccess.getElkPortAccess().getCommaKeyword_2_0()); }
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
	{ before(grammarAccess.getElkPortAccess().getPortElementParserRuleCall_2_1()); }
	rulePortElement
	{ after(grammarAccess.getElkPortAccess().getPortElementParserRuleCall_2_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__PortElement__Group_2__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__PortElement__Group_2__0__Impl
	rule__PortElement__Group_2__1
;
finally {
	restoreStackSize(stackSize);
}

rule__PortElement__Group_2__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getPortElementAccess().getKeyLabelsParserRuleCall_2_0()); }
	ruleKeyLabels
	{ after(grammarAccess.getPortElementAccess().getKeyLabelsParserRuleCall_2_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__PortElement__Group_2__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__PortElement__Group_2__1__Impl
	rule__PortElement__Group_2__2
;
finally {
	restoreStackSize(stackSize);
}

rule__PortElement__Group_2__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getPortElementAccess().getColonKeyword_2_1()); }
	':'
	{ after(grammarAccess.getPortElementAccess().getColonKeyword_2_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__PortElement__Group_2__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__PortElement__Group_2__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__PortElement__Group_2__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getPortElementAccess().getElkGraphElementLabelsParserRuleCall_2_2()); }
	ruleElkGraphElementLabels
	{ after(grammarAccess.getPortElementAccess().getElkGraphElementLabelsParserRuleCall_2_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__PortElement__Group_3__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__PortElement__Group_3__0__Impl
	rule__PortElement__Group_3__1
;
finally {
	restoreStackSize(stackSize);
}

rule__PortElement__Group_3__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getPortElementAccess().getKeyLayoutOptionsParserRuleCall_3_0()); }
	ruleKeyLayoutOptions
	{ after(grammarAccess.getPortElementAccess().getKeyLayoutOptionsParserRuleCall_3_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__PortElement__Group_3__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__PortElement__Group_3__1__Impl
	rule__PortElement__Group_3__2
;
finally {
	restoreStackSize(stackSize);
}

rule__PortElement__Group_3__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getPortElementAccess().getColonKeyword_3_1()); }
	':'
	{ after(grammarAccess.getPortElementAccess().getColonKeyword_3_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__PortElement__Group_3__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__PortElement__Group_3__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__PortElement__Group_3__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getPortElementAccess().getElkGraphElementPropertiesParserRuleCall_3_2()); }
	ruleElkGraphElementProperties
	{ after(grammarAccess.getPortElementAccess().getElkGraphElementPropertiesParserRuleCall_3_2()); }
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
	{ before(grammarAccess.getElkLabelAccess().getLeftCurlyBracketKeyword_0()); }
	'{'
	{ after(grammarAccess.getElkLabelAccess().getLeftCurlyBracketKeyword_0()); }
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
	{ before(grammarAccess.getElkLabelAccess().getLabelElementParserRuleCall_1()); }
	ruleLabelElement
	{ after(grammarAccess.getElkLabelAccess().getLabelElementParserRuleCall_1()); }
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
	{ before(grammarAccess.getElkLabelAccess().getGroup_2()); }
	(rule__ElkLabel__Group_2__0)*
	{ after(grammarAccess.getElkLabelAccess().getGroup_2()); }
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
	rule__ElkLabel__Group__4
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
	{ before(grammarAccess.getElkLabelAccess().getCommaKeyword_3()); }
	(',')?
	{ after(grammarAccess.getElkLabelAccess().getCommaKeyword_3()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkLabel__Group__4
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkLabel__Group__4__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkLabel__Group__4__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkLabelAccess().getRightCurlyBracketKeyword_4()); }
	'}'
	{ after(grammarAccess.getElkLabelAccess().getRightCurlyBracketKeyword_4()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkLabel__Group_2__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkLabel__Group_2__0__Impl
	rule__ElkLabel__Group_2__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkLabel__Group_2__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkLabelAccess().getCommaKeyword_2_0()); }
	','
	{ after(grammarAccess.getElkLabelAccess().getCommaKeyword_2_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkLabel__Group_2__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkLabel__Group_2__1__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkLabel__Group_2__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkLabelAccess().getLabelElementParserRuleCall_2_1()); }
	ruleLabelElement
	{ after(grammarAccess.getElkLabelAccess().getLabelElementParserRuleCall_2_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__LabelElement__Group_2__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__LabelElement__Group_2__0__Impl
	rule__LabelElement__Group_2__1
;
finally {
	restoreStackSize(stackSize);
}

rule__LabelElement__Group_2__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getLabelElementAccess().getKeyTextParserRuleCall_2_0()); }
	ruleKeyText
	{ after(grammarAccess.getLabelElementAccess().getKeyTextParserRuleCall_2_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__LabelElement__Group_2__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__LabelElement__Group_2__1__Impl
	rule__LabelElement__Group_2__2
;
finally {
	restoreStackSize(stackSize);
}

rule__LabelElement__Group_2__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getLabelElementAccess().getColonKeyword_2_1()); }
	':'
	{ after(grammarAccess.getLabelElementAccess().getColonKeyword_2_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__LabelElement__Group_2__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__LabelElement__Group_2__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__LabelElement__Group_2__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getLabelElementAccess().getTextAssignment_2_2()); }
	(rule__LabelElement__TextAssignment_2_2)
	{ after(grammarAccess.getLabelElementAccess().getTextAssignment_2_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__LabelElement__Group_3__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__LabelElement__Group_3__0__Impl
	rule__LabelElement__Group_3__1
;
finally {
	restoreStackSize(stackSize);
}

rule__LabelElement__Group_3__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getLabelElementAccess().getKeyLabelsParserRuleCall_3_0()); }
	ruleKeyLabels
	{ after(grammarAccess.getLabelElementAccess().getKeyLabelsParserRuleCall_3_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__LabelElement__Group_3__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__LabelElement__Group_3__1__Impl
	rule__LabelElement__Group_3__2
;
finally {
	restoreStackSize(stackSize);
}

rule__LabelElement__Group_3__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getLabelElementAccess().getColonKeyword_3_1()); }
	':'
	{ after(grammarAccess.getLabelElementAccess().getColonKeyword_3_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__LabelElement__Group_3__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__LabelElement__Group_3__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__LabelElement__Group_3__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getLabelElementAccess().getElkGraphElementLabelsParserRuleCall_3_2()); }
	ruleElkGraphElementLabels
	{ after(grammarAccess.getLabelElementAccess().getElkGraphElementLabelsParserRuleCall_3_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__LabelElement__Group_4__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__LabelElement__Group_4__0__Impl
	rule__LabelElement__Group_4__1
;
finally {
	restoreStackSize(stackSize);
}

rule__LabelElement__Group_4__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getLabelElementAccess().getKeyLayoutOptionsParserRuleCall_4_0()); }
	ruleKeyLayoutOptions
	{ after(grammarAccess.getLabelElementAccess().getKeyLayoutOptionsParserRuleCall_4_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__LabelElement__Group_4__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__LabelElement__Group_4__1__Impl
	rule__LabelElement__Group_4__2
;
finally {
	restoreStackSize(stackSize);
}

rule__LabelElement__Group_4__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getLabelElementAccess().getColonKeyword_4_1()); }
	':'
	{ after(grammarAccess.getLabelElementAccess().getColonKeyword_4_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__LabelElement__Group_4__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__LabelElement__Group_4__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__LabelElement__Group_4__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getLabelElementAccess().getElkGraphElementPropertiesParserRuleCall_4_2()); }
	ruleElkGraphElementProperties
	{ after(grammarAccess.getLabelElementAccess().getElkGraphElementPropertiesParserRuleCall_4_2()); }
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
	{ before(grammarAccess.getElkEdgeAccess().getLeftCurlyBracketKeyword_0()); }
	'{'
	{ after(grammarAccess.getElkEdgeAccess().getLeftCurlyBracketKeyword_0()); }
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
	{ before(grammarAccess.getElkEdgeAccess().getEdgeElementParserRuleCall_1()); }
	ruleEdgeElement
	{ after(grammarAccess.getElkEdgeAccess().getEdgeElementParserRuleCall_1()); }
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
	{ before(grammarAccess.getElkEdgeAccess().getGroup_2()); }
	(rule__ElkEdge__Group_2__0)*
	{ after(grammarAccess.getElkEdgeAccess().getGroup_2()); }
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
	{ before(grammarAccess.getElkEdgeAccess().getCommaKeyword_3()); }
	(',')?
	{ after(grammarAccess.getElkEdgeAccess().getCommaKeyword_3()); }
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
	{ before(grammarAccess.getElkEdgeAccess().getRightCurlyBracketKeyword_4()); }
	'}'
	{ after(grammarAccess.getElkEdgeAccess().getRightCurlyBracketKeyword_4()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkEdge__Group_2__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdge__Group_2__0__Impl
	rule__ElkEdge__Group_2__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdge__Group_2__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeAccess().getCommaKeyword_2_0()); }
	','
	{ after(grammarAccess.getElkEdgeAccess().getCommaKeyword_2_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdge__Group_2__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdge__Group_2__1__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdge__Group_2__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeAccess().getEdgeElementParserRuleCall_2_1()); }
	ruleEdgeElement
	{ after(grammarAccess.getElkEdgeAccess().getEdgeElementParserRuleCall_2_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__EdgeElement__Group_1__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__EdgeElement__Group_1__0__Impl
	rule__EdgeElement__Group_1__1
;
finally {
	restoreStackSize(stackSize);
}

rule__EdgeElement__Group_1__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getEdgeElementAccess().getKeySourcesParserRuleCall_1_0()); }
	ruleKeySources
	{ after(grammarAccess.getEdgeElementAccess().getKeySourcesParserRuleCall_1_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__EdgeElement__Group_1__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__EdgeElement__Group_1__1__Impl
	rule__EdgeElement__Group_1__2
;
finally {
	restoreStackSize(stackSize);
}

rule__EdgeElement__Group_1__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getEdgeElementAccess().getColonKeyword_1_1()); }
	':'
	{ after(grammarAccess.getEdgeElementAccess().getColonKeyword_1_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__EdgeElement__Group_1__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__EdgeElement__Group_1__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__EdgeElement__Group_1__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getEdgeElementAccess().getElkEdgeSourcesParserRuleCall_1_2()); }
	ruleElkEdgeSources
	{ after(grammarAccess.getEdgeElementAccess().getElkEdgeSourcesParserRuleCall_1_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__EdgeElement__Group_2__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__EdgeElement__Group_2__0__Impl
	rule__EdgeElement__Group_2__1
;
finally {
	restoreStackSize(stackSize);
}

rule__EdgeElement__Group_2__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getEdgeElementAccess().getKeyTargetsParserRuleCall_2_0()); }
	ruleKeyTargets
	{ after(grammarAccess.getEdgeElementAccess().getKeyTargetsParserRuleCall_2_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__EdgeElement__Group_2__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__EdgeElement__Group_2__1__Impl
	rule__EdgeElement__Group_2__2
;
finally {
	restoreStackSize(stackSize);
}

rule__EdgeElement__Group_2__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getEdgeElementAccess().getColonKeyword_2_1()); }
	':'
	{ after(grammarAccess.getEdgeElementAccess().getColonKeyword_2_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__EdgeElement__Group_2__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__EdgeElement__Group_2__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__EdgeElement__Group_2__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getEdgeElementAccess().getElkEdgeTargetsParserRuleCall_2_2()); }
	ruleElkEdgeTargets
	{ after(grammarAccess.getEdgeElementAccess().getElkEdgeTargetsParserRuleCall_2_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__EdgeElement__Group_3__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__EdgeElement__Group_3__0__Impl
	rule__EdgeElement__Group_3__1
;
finally {
	restoreStackSize(stackSize);
}

rule__EdgeElement__Group_3__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getEdgeElementAccess().getKeyLabelsParserRuleCall_3_0()); }
	ruleKeyLabels
	{ after(grammarAccess.getEdgeElementAccess().getKeyLabelsParserRuleCall_3_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__EdgeElement__Group_3__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__EdgeElement__Group_3__1__Impl
	rule__EdgeElement__Group_3__2
;
finally {
	restoreStackSize(stackSize);
}

rule__EdgeElement__Group_3__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getEdgeElementAccess().getColonKeyword_3_1()); }
	':'
	{ after(grammarAccess.getEdgeElementAccess().getColonKeyword_3_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__EdgeElement__Group_3__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__EdgeElement__Group_3__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__EdgeElement__Group_3__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getEdgeElementAccess().getElkGraphElementLabelsParserRuleCall_3_2()); }
	ruleElkGraphElementLabels
	{ after(grammarAccess.getEdgeElementAccess().getElkGraphElementLabelsParserRuleCall_3_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__EdgeElement__Group_4__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__EdgeElement__Group_4__0__Impl
	rule__EdgeElement__Group_4__1
;
finally {
	restoreStackSize(stackSize);
}

rule__EdgeElement__Group_4__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getEdgeElementAccess().getKeyLayoutOptionsParserRuleCall_4_0()); }
	ruleKeyLayoutOptions
	{ after(grammarAccess.getEdgeElementAccess().getKeyLayoutOptionsParserRuleCall_4_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__EdgeElement__Group_4__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__EdgeElement__Group_4__1__Impl
	rule__EdgeElement__Group_4__2
;
finally {
	restoreStackSize(stackSize);
}

rule__EdgeElement__Group_4__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getEdgeElementAccess().getColonKeyword_4_1()); }
	':'
	{ after(grammarAccess.getEdgeElementAccess().getColonKeyword_4_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__EdgeElement__Group_4__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__EdgeElement__Group_4__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__EdgeElement__Group_4__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getEdgeElementAccess().getElkGraphElementPropertiesParserRuleCall_4_2()); }
	ruleElkGraphElementProperties
	{ after(grammarAccess.getEdgeElementAccess().getElkGraphElementPropertiesParserRuleCall_4_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkEdgeSources__Group__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeSources__Group__0__Impl
	rule__ElkEdgeSources__Group__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSources__Group__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeSourcesAccess().getLeftSquareBracketKeyword_0()); }
	'['
	{ after(grammarAccess.getElkEdgeSourcesAccess().getLeftSquareBracketKeyword_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSources__Group__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeSources__Group__1__Impl
	rule__ElkEdgeSources__Group__2
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSources__Group__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeSourcesAccess().getGroup_1()); }
	(rule__ElkEdgeSources__Group_1__0)?
	{ after(grammarAccess.getElkEdgeSourcesAccess().getGroup_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSources__Group__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeSources__Group__2__Impl
	rule__ElkEdgeSources__Group__3
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSources__Group__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeSourcesAccess().getCommaKeyword_2()); }
	(',')?
	{ after(grammarAccess.getElkEdgeSourcesAccess().getCommaKeyword_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSources__Group__3
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeSources__Group__3__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSources__Group__3__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeSourcesAccess().getRightSquareBracketKeyword_3()); }
	']'
	{ after(grammarAccess.getElkEdgeSourcesAccess().getRightSquareBracketKeyword_3()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkEdgeSources__Group_1__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeSources__Group_1__0__Impl
	rule__ElkEdgeSources__Group_1__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSources__Group_1__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeSourcesAccess().getSourcesAssignment_1_0()); }
	(rule__ElkEdgeSources__SourcesAssignment_1_0)
	{ after(grammarAccess.getElkEdgeSourcesAccess().getSourcesAssignment_1_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSources__Group_1__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeSources__Group_1__1__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSources__Group_1__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeSourcesAccess().getGroup_1_1()); }
	(rule__ElkEdgeSources__Group_1_1__0)*
	{ after(grammarAccess.getElkEdgeSourcesAccess().getGroup_1_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkEdgeSources__Group_1_1__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeSources__Group_1_1__0__Impl
	rule__ElkEdgeSources__Group_1_1__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSources__Group_1_1__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeSourcesAccess().getCommaKeyword_1_1_0()); }
	','
	{ after(grammarAccess.getElkEdgeSourcesAccess().getCommaKeyword_1_1_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSources__Group_1_1__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeSources__Group_1_1__1__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSources__Group_1_1__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeSourcesAccess().getSourcesAssignment_1_1_1()); }
	(rule__ElkEdgeSources__SourcesAssignment_1_1_1)
	{ after(grammarAccess.getElkEdgeSourcesAccess().getSourcesAssignment_1_1_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkEdgeTargets__Group__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeTargets__Group__0__Impl
	rule__ElkEdgeTargets__Group__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeTargets__Group__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeTargetsAccess().getLeftSquareBracketKeyword_0()); }
	'['
	{ after(grammarAccess.getElkEdgeTargetsAccess().getLeftSquareBracketKeyword_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeTargets__Group__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeTargets__Group__1__Impl
	rule__ElkEdgeTargets__Group__2
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeTargets__Group__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeTargetsAccess().getGroup_1()); }
	(rule__ElkEdgeTargets__Group_1__0)?
	{ after(grammarAccess.getElkEdgeTargetsAccess().getGroup_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeTargets__Group__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeTargets__Group__2__Impl
	rule__ElkEdgeTargets__Group__3
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeTargets__Group__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeTargetsAccess().getCommaKeyword_2()); }
	(',')?
	{ after(grammarAccess.getElkEdgeTargetsAccess().getCommaKeyword_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeTargets__Group__3
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeTargets__Group__3__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeTargets__Group__3__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeTargetsAccess().getRightSquareBracketKeyword_3()); }
	']'
	{ after(grammarAccess.getElkEdgeTargetsAccess().getRightSquareBracketKeyword_3()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkEdgeTargets__Group_1__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeTargets__Group_1__0__Impl
	rule__ElkEdgeTargets__Group_1__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeTargets__Group_1__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeTargetsAccess().getTargetsAssignment_1_0()); }
	(rule__ElkEdgeTargets__TargetsAssignment_1_0)
	{ after(grammarAccess.getElkEdgeTargetsAccess().getTargetsAssignment_1_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeTargets__Group_1__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeTargets__Group_1__1__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeTargets__Group_1__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeTargetsAccess().getGroup_1_1()); }
	(rule__ElkEdgeTargets__Group_1_1__0)*
	{ after(grammarAccess.getElkEdgeTargetsAccess().getGroup_1_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkEdgeTargets__Group_1_1__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeTargets__Group_1_1__0__Impl
	rule__ElkEdgeTargets__Group_1_1__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeTargets__Group_1_1__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeTargetsAccess().getCommaKeyword_1_1_0()); }
	','
	{ after(grammarAccess.getElkEdgeTargetsAccess().getCommaKeyword_1_1_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeTargets__Group_1_1__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkEdgeTargets__Group_1_1__1__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeTargets__Group_1_1__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkEdgeTargetsAccess().getTargetsAssignment_1_1_1()); }
	(rule__ElkEdgeTargets__TargetsAssignment_1_1_1)
	{ after(grammarAccess.getElkEdgeTargetsAccess().getTargetsAssignment_1_1_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkId__Group__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkId__Group__0__Impl
	rule__ElkId__Group__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkId__Group__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkIdAccess().getKeyIdParserRuleCall_0()); }
	ruleKeyId
	{ after(grammarAccess.getElkIdAccess().getKeyIdParserRuleCall_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkId__Group__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkId__Group__1__Impl
	rule__ElkId__Group__2
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkId__Group__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkIdAccess().getColonKeyword_1()); }
	':'
	{ after(grammarAccess.getElkIdAccess().getColonKeyword_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkId__Group__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkId__Group__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkId__Group__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkIdAccess().getIdentifierAssignment_2()); }
	(rule__ElkId__IdentifierAssignment_2)
	{ after(grammarAccess.getElkIdAccess().getIdentifierAssignment_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkNodeChildren__Group__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkNodeChildren__Group__0__Impl
	rule__ElkNodeChildren__Group__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNodeChildren__Group__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkNodeChildrenAccess().getLeftSquareBracketKeyword_0()); }
	'['
	{ after(grammarAccess.getElkNodeChildrenAccess().getLeftSquareBracketKeyword_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNodeChildren__Group__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkNodeChildren__Group__1__Impl
	rule__ElkNodeChildren__Group__2
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNodeChildren__Group__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkNodeChildrenAccess().getGroup_1()); }
	(rule__ElkNodeChildren__Group_1__0)?
	{ after(grammarAccess.getElkNodeChildrenAccess().getGroup_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNodeChildren__Group__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkNodeChildren__Group__2__Impl
	rule__ElkNodeChildren__Group__3
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNodeChildren__Group__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkNodeChildrenAccess().getCommaKeyword_2()); }
	(',')?
	{ after(grammarAccess.getElkNodeChildrenAccess().getCommaKeyword_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNodeChildren__Group__3
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkNodeChildren__Group__3__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNodeChildren__Group__3__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkNodeChildrenAccess().getRightSquareBracketKeyword_3()); }
	']'
	{ after(grammarAccess.getElkNodeChildrenAccess().getRightSquareBracketKeyword_3()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkNodeChildren__Group_1__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkNodeChildren__Group_1__0__Impl
	rule__ElkNodeChildren__Group_1__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNodeChildren__Group_1__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkNodeChildrenAccess().getChildrenAssignment_1_0()); }
	(rule__ElkNodeChildren__ChildrenAssignment_1_0)
	{ after(grammarAccess.getElkNodeChildrenAccess().getChildrenAssignment_1_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNodeChildren__Group_1__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkNodeChildren__Group_1__1__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNodeChildren__Group_1__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkNodeChildrenAccess().getGroup_1_1()); }
	(rule__ElkNodeChildren__Group_1_1__0)*
	{ after(grammarAccess.getElkNodeChildrenAccess().getGroup_1_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkNodeChildren__Group_1_1__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkNodeChildren__Group_1_1__0__Impl
	rule__ElkNodeChildren__Group_1_1__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNodeChildren__Group_1_1__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkNodeChildrenAccess().getCommaKeyword_1_1_0()); }
	','
	{ after(grammarAccess.getElkNodeChildrenAccess().getCommaKeyword_1_1_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNodeChildren__Group_1_1__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkNodeChildren__Group_1_1__1__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNodeChildren__Group_1_1__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkNodeChildrenAccess().getChildrenAssignment_1_1_1()); }
	(rule__ElkNodeChildren__ChildrenAssignment_1_1_1)
	{ after(grammarAccess.getElkNodeChildrenAccess().getChildrenAssignment_1_1_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkNodePorts__Group__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkNodePorts__Group__0__Impl
	rule__ElkNodePorts__Group__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNodePorts__Group__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkNodePortsAccess().getLeftSquareBracketKeyword_0()); }
	'['
	{ after(grammarAccess.getElkNodePortsAccess().getLeftSquareBracketKeyword_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNodePorts__Group__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkNodePorts__Group__1__Impl
	rule__ElkNodePorts__Group__2
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNodePorts__Group__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkNodePortsAccess().getGroup_1()); }
	(rule__ElkNodePorts__Group_1__0)?
	{ after(grammarAccess.getElkNodePortsAccess().getGroup_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNodePorts__Group__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkNodePorts__Group__2__Impl
	rule__ElkNodePorts__Group__3
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNodePorts__Group__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkNodePortsAccess().getCommaKeyword_2()); }
	(',')?
	{ after(grammarAccess.getElkNodePortsAccess().getCommaKeyword_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNodePorts__Group__3
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkNodePorts__Group__3__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNodePorts__Group__3__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkNodePortsAccess().getRightSquareBracketKeyword_3()); }
	']'
	{ after(grammarAccess.getElkNodePortsAccess().getRightSquareBracketKeyword_3()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkNodePorts__Group_1__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkNodePorts__Group_1__0__Impl
	rule__ElkNodePorts__Group_1__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNodePorts__Group_1__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkNodePortsAccess().getPortsAssignment_1_0()); }
	(rule__ElkNodePorts__PortsAssignment_1_0)
	{ after(grammarAccess.getElkNodePortsAccess().getPortsAssignment_1_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNodePorts__Group_1__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkNodePorts__Group_1__1__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNodePorts__Group_1__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkNodePortsAccess().getGroup_1_1()); }
	(rule__ElkNodePorts__Group_1_1__0)*
	{ after(grammarAccess.getElkNodePortsAccess().getGroup_1_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkNodePorts__Group_1_1__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkNodePorts__Group_1_1__0__Impl
	rule__ElkNodePorts__Group_1_1__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNodePorts__Group_1_1__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkNodePortsAccess().getCommaKeyword_1_1_0()); }
	','
	{ after(grammarAccess.getElkNodePortsAccess().getCommaKeyword_1_1_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNodePorts__Group_1_1__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkNodePorts__Group_1_1__1__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNodePorts__Group_1_1__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkNodePortsAccess().getPortsAssignment_1_1_1()); }
	(rule__ElkNodePorts__PortsAssignment_1_1_1)
	{ after(grammarAccess.getElkNodePortsAccess().getPortsAssignment_1_1_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkNodeEdges__Group__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkNodeEdges__Group__0__Impl
	rule__ElkNodeEdges__Group__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNodeEdges__Group__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkNodeEdgesAccess().getLeftSquareBracketKeyword_0()); }
	'['
	{ after(grammarAccess.getElkNodeEdgesAccess().getLeftSquareBracketKeyword_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNodeEdges__Group__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkNodeEdges__Group__1__Impl
	rule__ElkNodeEdges__Group__2
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNodeEdges__Group__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkNodeEdgesAccess().getGroup_1()); }
	(rule__ElkNodeEdges__Group_1__0)?
	{ after(grammarAccess.getElkNodeEdgesAccess().getGroup_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNodeEdges__Group__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkNodeEdges__Group__2__Impl
	rule__ElkNodeEdges__Group__3
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNodeEdges__Group__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkNodeEdgesAccess().getCommaKeyword_2()); }
	(',')?
	{ after(grammarAccess.getElkNodeEdgesAccess().getCommaKeyword_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNodeEdges__Group__3
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkNodeEdges__Group__3__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNodeEdges__Group__3__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkNodeEdgesAccess().getRightSquareBracketKeyword_3()); }
	']'
	{ after(grammarAccess.getElkNodeEdgesAccess().getRightSquareBracketKeyword_3()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkNodeEdges__Group_1__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkNodeEdges__Group_1__0__Impl
	rule__ElkNodeEdges__Group_1__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNodeEdges__Group_1__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkNodeEdgesAccess().getContainedEdgesAssignment_1_0()); }
	(rule__ElkNodeEdges__ContainedEdgesAssignment_1_0)
	{ after(grammarAccess.getElkNodeEdgesAccess().getContainedEdgesAssignment_1_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNodeEdges__Group_1__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkNodeEdges__Group_1__1__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNodeEdges__Group_1__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkNodeEdgesAccess().getGroup_1_1()); }
	(rule__ElkNodeEdges__Group_1_1__0)*
	{ after(grammarAccess.getElkNodeEdgesAccess().getGroup_1_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkNodeEdges__Group_1_1__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkNodeEdges__Group_1_1__0__Impl
	rule__ElkNodeEdges__Group_1_1__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNodeEdges__Group_1_1__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkNodeEdgesAccess().getCommaKeyword_1_1_0()); }
	','
	{ after(grammarAccess.getElkNodeEdgesAccess().getCommaKeyword_1_1_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNodeEdges__Group_1_1__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkNodeEdges__Group_1_1__1__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNodeEdges__Group_1_1__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkNodeEdgesAccess().getContainedEdgesAssignment_1_1_1()); }
	(rule__ElkNodeEdges__ContainedEdgesAssignment_1_1_1)
	{ after(grammarAccess.getElkNodeEdgesAccess().getContainedEdgesAssignment_1_1_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkGraphElementLabels__Group__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkGraphElementLabels__Group__0__Impl
	rule__ElkGraphElementLabels__Group__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkGraphElementLabels__Group__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkGraphElementLabelsAccess().getLeftSquareBracketKeyword_0()); }
	'['
	{ after(grammarAccess.getElkGraphElementLabelsAccess().getLeftSquareBracketKeyword_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkGraphElementLabels__Group__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkGraphElementLabels__Group__1__Impl
	rule__ElkGraphElementLabels__Group__2
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkGraphElementLabels__Group__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkGraphElementLabelsAccess().getGroup_1()); }
	(rule__ElkGraphElementLabels__Group_1__0)?
	{ after(grammarAccess.getElkGraphElementLabelsAccess().getGroup_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkGraphElementLabels__Group__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkGraphElementLabels__Group__2__Impl
	rule__ElkGraphElementLabels__Group__3
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkGraphElementLabels__Group__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkGraphElementLabelsAccess().getCommaKeyword_2()); }
	(',')?
	{ after(grammarAccess.getElkGraphElementLabelsAccess().getCommaKeyword_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkGraphElementLabels__Group__3
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkGraphElementLabels__Group__3__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkGraphElementLabels__Group__3__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkGraphElementLabelsAccess().getRightSquareBracketKeyword_3()); }
	']'
	{ after(grammarAccess.getElkGraphElementLabelsAccess().getRightSquareBracketKeyword_3()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkGraphElementLabels__Group_1__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkGraphElementLabels__Group_1__0__Impl
	rule__ElkGraphElementLabels__Group_1__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkGraphElementLabels__Group_1__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkGraphElementLabelsAccess().getLabelsAssignment_1_0()); }
	(rule__ElkGraphElementLabels__LabelsAssignment_1_0)
	{ after(grammarAccess.getElkGraphElementLabelsAccess().getLabelsAssignment_1_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkGraphElementLabels__Group_1__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkGraphElementLabels__Group_1__1__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkGraphElementLabels__Group_1__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkGraphElementLabelsAccess().getGroup_1_1()); }
	(rule__ElkGraphElementLabels__Group_1_1__0)*
	{ after(grammarAccess.getElkGraphElementLabelsAccess().getGroup_1_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkGraphElementLabels__Group_1_1__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkGraphElementLabels__Group_1_1__0__Impl
	rule__ElkGraphElementLabels__Group_1_1__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkGraphElementLabels__Group_1_1__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkGraphElementLabelsAccess().getCommaKeyword_1_1_0()); }
	','
	{ after(grammarAccess.getElkGraphElementLabelsAccess().getCommaKeyword_1_1_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkGraphElementLabels__Group_1_1__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkGraphElementLabels__Group_1_1__1__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkGraphElementLabels__Group_1_1__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkGraphElementLabelsAccess().getLabelsAssignment_1_1_1()); }
	(rule__ElkGraphElementLabels__LabelsAssignment_1_1_1)
	{ after(grammarAccess.getElkGraphElementLabelsAccess().getLabelsAssignment_1_1_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkGraphElementProperties__Group__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkGraphElementProperties__Group__0__Impl
	rule__ElkGraphElementProperties__Group__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkGraphElementProperties__Group__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkGraphElementPropertiesAccess().getLeftCurlyBracketKeyword_0()); }
	'{'
	{ after(grammarAccess.getElkGraphElementPropertiesAccess().getLeftCurlyBracketKeyword_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkGraphElementProperties__Group__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkGraphElementProperties__Group__1__Impl
	rule__ElkGraphElementProperties__Group__2
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkGraphElementProperties__Group__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkGraphElementPropertiesAccess().getGroup_1()); }
	(rule__ElkGraphElementProperties__Group_1__0)?
	{ after(grammarAccess.getElkGraphElementPropertiesAccess().getGroup_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkGraphElementProperties__Group__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkGraphElementProperties__Group__2__Impl
	rule__ElkGraphElementProperties__Group__3
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkGraphElementProperties__Group__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkGraphElementPropertiesAccess().getCommaKeyword_2()); }
	(',')?
	{ after(grammarAccess.getElkGraphElementPropertiesAccess().getCommaKeyword_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkGraphElementProperties__Group__3
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkGraphElementProperties__Group__3__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkGraphElementProperties__Group__3__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkGraphElementPropertiesAccess().getRightCurlyBracketKeyword_3()); }
	'}'
	{ after(grammarAccess.getElkGraphElementPropertiesAccess().getRightCurlyBracketKeyword_3()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkGraphElementProperties__Group_1__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkGraphElementProperties__Group_1__0__Impl
	rule__ElkGraphElementProperties__Group_1__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkGraphElementProperties__Group_1__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkGraphElementPropertiesAccess().getPropertiesAssignment_1_0()); }
	(rule__ElkGraphElementProperties__PropertiesAssignment_1_0)
	{ after(grammarAccess.getElkGraphElementPropertiesAccess().getPropertiesAssignment_1_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkGraphElementProperties__Group_1__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkGraphElementProperties__Group_1__1__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkGraphElementProperties__Group_1__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkGraphElementPropertiesAccess().getGroup_1_1()); }
	(rule__ElkGraphElementProperties__Group_1_1__0)*
	{ after(grammarAccess.getElkGraphElementPropertiesAccess().getGroup_1_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ElkGraphElementProperties__Group_1_1__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkGraphElementProperties__Group_1_1__0__Impl
	rule__ElkGraphElementProperties__Group_1_1__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkGraphElementProperties__Group_1_1__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkGraphElementPropertiesAccess().getCommaKeyword_1_1_0()); }
	','
	{ after(grammarAccess.getElkGraphElementPropertiesAccess().getCommaKeyword_1_1_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkGraphElementProperties__Group_1_1__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ElkGraphElementProperties__Group_1_1__1__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkGraphElementProperties__Group_1_1__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getElkGraphElementPropertiesAccess().getPropertiesAssignment_1_1_1()); }
	(rule__ElkGraphElementProperties__PropertiesAssignment_1_1_1)
	{ after(grammarAccess.getElkGraphElementPropertiesAccess().getPropertiesAssignment_1_1_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ShapeElement__Group_0__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ShapeElement__Group_0__0__Impl
	rule__ShapeElement__Group_0__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeElement__Group_0__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getShapeElementAccess().getKeyXParserRuleCall_0_0()); }
	ruleKeyX
	{ after(grammarAccess.getShapeElementAccess().getKeyXParserRuleCall_0_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeElement__Group_0__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ShapeElement__Group_0__1__Impl
	rule__ShapeElement__Group_0__2
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeElement__Group_0__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getShapeElementAccess().getColonKeyword_0_1()); }
	':'
	{ after(grammarAccess.getShapeElementAccess().getColonKeyword_0_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeElement__Group_0__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ShapeElement__Group_0__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeElement__Group_0__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getShapeElementAccess().getXAssignment_0_2()); }
	(rule__ShapeElement__XAssignment_0_2)
	{ after(grammarAccess.getShapeElementAccess().getXAssignment_0_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ShapeElement__Group_1__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ShapeElement__Group_1__0__Impl
	rule__ShapeElement__Group_1__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeElement__Group_1__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getShapeElementAccess().getKeyYParserRuleCall_1_0()); }
	ruleKeyY
	{ after(grammarAccess.getShapeElementAccess().getKeyYParserRuleCall_1_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeElement__Group_1__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ShapeElement__Group_1__1__Impl
	rule__ShapeElement__Group_1__2
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeElement__Group_1__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getShapeElementAccess().getColonKeyword_1_1()); }
	':'
	{ after(grammarAccess.getShapeElementAccess().getColonKeyword_1_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeElement__Group_1__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ShapeElement__Group_1__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeElement__Group_1__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getShapeElementAccess().getYAssignment_1_2()); }
	(rule__ShapeElement__YAssignment_1_2)
	{ after(grammarAccess.getShapeElementAccess().getYAssignment_1_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ShapeElement__Group_2__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ShapeElement__Group_2__0__Impl
	rule__ShapeElement__Group_2__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeElement__Group_2__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getShapeElementAccess().getKeyWidthParserRuleCall_2_0()); }
	ruleKeyWidth
	{ after(grammarAccess.getShapeElementAccess().getKeyWidthParserRuleCall_2_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeElement__Group_2__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ShapeElement__Group_2__1__Impl
	rule__ShapeElement__Group_2__2
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeElement__Group_2__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getShapeElementAccess().getColonKeyword_2_1()); }
	':'
	{ after(grammarAccess.getShapeElementAccess().getColonKeyword_2_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeElement__Group_2__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ShapeElement__Group_2__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeElement__Group_2__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getShapeElementAccess().getWidthAssignment_2_2()); }
	(rule__ShapeElement__WidthAssignment_2_2)
	{ after(grammarAccess.getShapeElementAccess().getWidthAssignment_2_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__ShapeElement__Group_3__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ShapeElement__Group_3__0__Impl
	rule__ShapeElement__Group_3__1
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeElement__Group_3__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getShapeElementAccess().getKeyHeightParserRuleCall_3_0()); }
	ruleKeyHeight
	{ after(grammarAccess.getShapeElementAccess().getKeyHeightParserRuleCall_3_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeElement__Group_3__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ShapeElement__Group_3__1__Impl
	rule__ShapeElement__Group_3__2
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeElement__Group_3__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getShapeElementAccess().getColonKeyword_3_1()); }
	':'
	{ after(grammarAccess.getShapeElementAccess().getColonKeyword_3_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeElement__Group_3__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__ShapeElement__Group_3__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeElement__Group_3__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getShapeElementAccess().getHeightAssignment_3_2()); }
	(rule__ShapeElement__HeightAssignment_3_2)
	{ after(grammarAccess.getShapeElementAccess().getHeightAssignment_3_2()); }
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


rule__JsonObject__Group__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__JsonObject__Group__0__Impl
	rule__JsonObject__Group__1
;
finally {
	restoreStackSize(stackSize);
}

rule__JsonObject__Group__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getJsonObjectAccess().getLeftCurlyBracketKeyword_0()); }
	'{'
	{ after(grammarAccess.getJsonObjectAccess().getLeftCurlyBracketKeyword_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__JsonObject__Group__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__JsonObject__Group__1__Impl
	rule__JsonObject__Group__2
;
finally {
	restoreStackSize(stackSize);
}

rule__JsonObject__Group__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getJsonObjectAccess().getGroup_1()); }
	(rule__JsonObject__Group_1__0)?
	{ after(grammarAccess.getJsonObjectAccess().getGroup_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__JsonObject__Group__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__JsonObject__Group__2__Impl
	rule__JsonObject__Group__3
;
finally {
	restoreStackSize(stackSize);
}

rule__JsonObject__Group__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getJsonObjectAccess().getCommaKeyword_2()); }
	(',')?
	{ after(grammarAccess.getJsonObjectAccess().getCommaKeyword_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__JsonObject__Group__3
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__JsonObject__Group__3__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__JsonObject__Group__3__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getJsonObjectAccess().getRightCurlyBracketKeyword_3()); }
	'}'
	{ after(grammarAccess.getJsonObjectAccess().getRightCurlyBracketKeyword_3()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__JsonObject__Group_1__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__JsonObject__Group_1__0__Impl
	rule__JsonObject__Group_1__1
;
finally {
	restoreStackSize(stackSize);
}

rule__JsonObject__Group_1__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getJsonObjectAccess().getJsonMemberParserRuleCall_1_0()); }
	ruleJsonMember
	{ after(grammarAccess.getJsonObjectAccess().getJsonMemberParserRuleCall_1_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__JsonObject__Group_1__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__JsonObject__Group_1__1__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__JsonObject__Group_1__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getJsonObjectAccess().getGroup_1_1()); }
	(rule__JsonObject__Group_1_1__0)*
	{ after(grammarAccess.getJsonObjectAccess().getGroup_1_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__JsonObject__Group_1_1__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__JsonObject__Group_1_1__0__Impl
	rule__JsonObject__Group_1_1__1
;
finally {
	restoreStackSize(stackSize);
}

rule__JsonObject__Group_1_1__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getJsonObjectAccess().getCommaKeyword_1_1_0()); }
	','
	{ after(grammarAccess.getJsonObjectAccess().getCommaKeyword_1_1_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__JsonObject__Group_1_1__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__JsonObject__Group_1_1__1__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__JsonObject__Group_1_1__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getJsonObjectAccess().getJsonMemberParserRuleCall_1_1_1()); }
	ruleJsonMember
	{ after(grammarAccess.getJsonObjectAccess().getJsonMemberParserRuleCall_1_1_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__JsonArray__Group__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__JsonArray__Group__0__Impl
	rule__JsonArray__Group__1
;
finally {
	restoreStackSize(stackSize);
}

rule__JsonArray__Group__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getJsonArrayAccess().getLeftSquareBracketKeyword_0()); }
	'['
	{ after(grammarAccess.getJsonArrayAccess().getLeftSquareBracketKeyword_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__JsonArray__Group__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__JsonArray__Group__1__Impl
	rule__JsonArray__Group__2
;
finally {
	restoreStackSize(stackSize);
}

rule__JsonArray__Group__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getJsonArrayAccess().getGroup_1()); }
	(rule__JsonArray__Group_1__0)?
	{ after(grammarAccess.getJsonArrayAccess().getGroup_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__JsonArray__Group__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__JsonArray__Group__2__Impl
	rule__JsonArray__Group__3
;
finally {
	restoreStackSize(stackSize);
}

rule__JsonArray__Group__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getJsonArrayAccess().getCommaKeyword_2()); }
	(',')?
	{ after(grammarAccess.getJsonArrayAccess().getCommaKeyword_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__JsonArray__Group__3
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__JsonArray__Group__3__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__JsonArray__Group__3__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getJsonArrayAccess().getRightSquareBracketKeyword_3()); }
	']'
	{ after(grammarAccess.getJsonArrayAccess().getRightSquareBracketKeyword_3()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__JsonArray__Group_1__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__JsonArray__Group_1__0__Impl
	rule__JsonArray__Group_1__1
;
finally {
	restoreStackSize(stackSize);
}

rule__JsonArray__Group_1__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getJsonArrayAccess().getJsonValueParserRuleCall_1_0()); }
	ruleJsonValue
	{ after(grammarAccess.getJsonArrayAccess().getJsonValueParserRuleCall_1_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__JsonArray__Group_1__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__JsonArray__Group_1__1__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__JsonArray__Group_1__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getJsonArrayAccess().getGroup_1_1()); }
	(rule__JsonArray__Group_1_1__0)*
	{ after(grammarAccess.getJsonArrayAccess().getGroup_1_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__JsonArray__Group_1_1__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__JsonArray__Group_1_1__0__Impl
	rule__JsonArray__Group_1_1__1
;
finally {
	restoreStackSize(stackSize);
}

rule__JsonArray__Group_1_1__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getJsonArrayAccess().getCommaKeyword_1_1_0()); }
	','
	{ after(grammarAccess.getJsonArrayAccess().getCommaKeyword_1_1_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__JsonArray__Group_1_1__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__JsonArray__Group_1_1__1__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__JsonArray__Group_1_1__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getJsonArrayAccess().getJsonValueParserRuleCall_1_1_1()); }
	ruleJsonValue
	{ after(grammarAccess.getJsonArrayAccess().getJsonValueParserRuleCall_1_1_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__JsonMember__Group__0
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__JsonMember__Group__0__Impl
	rule__JsonMember__Group__1
;
finally {
	restoreStackSize(stackSize);
}

rule__JsonMember__Group__0__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getJsonMemberAccess().getAlternatives_0()); }
	(rule__JsonMember__Alternatives_0)
	{ after(grammarAccess.getJsonMemberAccess().getAlternatives_0()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__JsonMember__Group__1
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__JsonMember__Group__1__Impl
	rule__JsonMember__Group__2
;
finally {
	restoreStackSize(stackSize);
}

rule__JsonMember__Group__1__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getJsonMemberAccess().getColonKeyword_1()); }
	':'
	{ after(grammarAccess.getJsonMemberAccess().getColonKeyword_1()); }
)
;
finally {
	restoreStackSize(stackSize);
}

rule__JsonMember__Group__2
	@init {
		int stackSize = keepStackSize();
	}
:
	rule__JsonMember__Group__2__Impl
;
finally {
	restoreStackSize(stackSize);
}

rule__JsonMember__Group__2__Impl
	@init {
		int stackSize = keepStackSize();
	}
:
(
	{ before(grammarAccess.getJsonMemberAccess().getJsonValueParserRuleCall_2()); }
	ruleJsonValue
	{ after(grammarAccess.getJsonMemberAccess().getJsonValueParserRuleCall_2()); }
)
;
finally {
	restoreStackSize(stackSize);
}


rule__LabelElement__TextAssignment_2_2
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getLabelElementAccess().getTextSTRINGTerminalRuleCall_2_2_0()); }
		RULE_STRING
		{ after(grammarAccess.getLabelElementAccess().getTextSTRINGTerminalRuleCall_2_2_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSources__SourcesAssignment_1_0
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkEdgeSourcesAccess().getSourcesElkConnectableShapeCrossReference_1_0_0()); }
		(
			{ before(grammarAccess.getElkEdgeSourcesAccess().getSourcesElkConnectableShapeSTRINGTerminalRuleCall_1_0_0_1()); }
			RULE_STRING
			{ after(grammarAccess.getElkEdgeSourcesAccess().getSourcesElkConnectableShapeSTRINGTerminalRuleCall_1_0_0_1()); }
		)
		{ after(grammarAccess.getElkEdgeSourcesAccess().getSourcesElkConnectableShapeCrossReference_1_0_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeSources__SourcesAssignment_1_1_1
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkEdgeSourcesAccess().getSourcesElkConnectableShapeCrossReference_1_1_1_0()); }
		(
			{ before(grammarAccess.getElkEdgeSourcesAccess().getSourcesElkConnectableShapeSTRINGTerminalRuleCall_1_1_1_0_1()); }
			RULE_STRING
			{ after(grammarAccess.getElkEdgeSourcesAccess().getSourcesElkConnectableShapeSTRINGTerminalRuleCall_1_1_1_0_1()); }
		)
		{ after(grammarAccess.getElkEdgeSourcesAccess().getSourcesElkConnectableShapeCrossReference_1_1_1_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeTargets__TargetsAssignment_1_0
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkEdgeTargetsAccess().getTargetsElkConnectableShapeCrossReference_1_0_0()); }
		(
			{ before(grammarAccess.getElkEdgeTargetsAccess().getTargetsElkConnectableShapeSTRINGTerminalRuleCall_1_0_0_1()); }
			RULE_STRING
			{ after(grammarAccess.getElkEdgeTargetsAccess().getTargetsElkConnectableShapeSTRINGTerminalRuleCall_1_0_0_1()); }
		)
		{ after(grammarAccess.getElkEdgeTargetsAccess().getTargetsElkConnectableShapeCrossReference_1_0_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkEdgeTargets__TargetsAssignment_1_1_1
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkEdgeTargetsAccess().getTargetsElkConnectableShapeCrossReference_1_1_1_0()); }
		(
			{ before(grammarAccess.getElkEdgeTargetsAccess().getTargetsElkConnectableShapeSTRINGTerminalRuleCall_1_1_1_0_1()); }
			RULE_STRING
			{ after(grammarAccess.getElkEdgeTargetsAccess().getTargetsElkConnectableShapeSTRINGTerminalRuleCall_1_1_1_0_1()); }
		)
		{ after(grammarAccess.getElkEdgeTargetsAccess().getTargetsElkConnectableShapeCrossReference_1_1_1_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkId__IdentifierAssignment_2
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkIdAccess().getIdentifierSTRINGTerminalRuleCall_2_0()); }
		RULE_STRING
		{ after(grammarAccess.getElkIdAccess().getIdentifierSTRINGTerminalRuleCall_2_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNodeChildren__ChildrenAssignment_1_0
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkNodeChildrenAccess().getChildrenElkNodeParserRuleCall_1_0_0()); }
		ruleElkNode
		{ after(grammarAccess.getElkNodeChildrenAccess().getChildrenElkNodeParserRuleCall_1_0_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNodeChildren__ChildrenAssignment_1_1_1
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkNodeChildrenAccess().getChildrenElkNodeParserRuleCall_1_1_1_0()); }
		ruleElkNode
		{ after(grammarAccess.getElkNodeChildrenAccess().getChildrenElkNodeParserRuleCall_1_1_1_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNodePorts__PortsAssignment_1_0
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkNodePortsAccess().getPortsElkPortParserRuleCall_1_0_0()); }
		ruleElkPort
		{ after(grammarAccess.getElkNodePortsAccess().getPortsElkPortParserRuleCall_1_0_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNodePorts__PortsAssignment_1_1_1
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkNodePortsAccess().getPortsElkPortParserRuleCall_1_1_1_0()); }
		ruleElkPort
		{ after(grammarAccess.getElkNodePortsAccess().getPortsElkPortParserRuleCall_1_1_1_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNodeEdges__ContainedEdgesAssignment_1_0
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkNodeEdgesAccess().getContainedEdgesElkEdgeParserRuleCall_1_0_0()); }
		ruleElkEdge
		{ after(grammarAccess.getElkNodeEdgesAccess().getContainedEdgesElkEdgeParserRuleCall_1_0_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkNodeEdges__ContainedEdgesAssignment_1_1_1
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkNodeEdgesAccess().getContainedEdgesElkEdgeParserRuleCall_1_1_1_0()); }
		ruleElkEdge
		{ after(grammarAccess.getElkNodeEdgesAccess().getContainedEdgesElkEdgeParserRuleCall_1_1_1_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkGraphElementLabels__LabelsAssignment_1_0
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkGraphElementLabelsAccess().getLabelsElkLabelParserRuleCall_1_0_0()); }
		ruleElkLabel
		{ after(grammarAccess.getElkGraphElementLabelsAccess().getLabelsElkLabelParserRuleCall_1_0_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkGraphElementLabels__LabelsAssignment_1_1_1
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkGraphElementLabelsAccess().getLabelsElkLabelParserRuleCall_1_1_1_0()); }
		ruleElkLabel
		{ after(grammarAccess.getElkGraphElementLabelsAccess().getLabelsElkLabelParserRuleCall_1_1_1_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkGraphElementProperties__PropertiesAssignment_1_0
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkGraphElementPropertiesAccess().getPropertiesPropertyParserRuleCall_1_0_0()); }
		ruleProperty
		{ after(grammarAccess.getElkGraphElementPropertiesAccess().getPropertiesPropertyParserRuleCall_1_0_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ElkGraphElementProperties__PropertiesAssignment_1_1_1
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getElkGraphElementPropertiesAccess().getPropertiesPropertyParserRuleCall_1_1_1_0()); }
		ruleProperty
		{ after(grammarAccess.getElkGraphElementPropertiesAccess().getPropertiesPropertyParserRuleCall_1_1_1_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeElement__XAssignment_0_2
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getShapeElementAccess().getXNumberParserRuleCall_0_2_0()); }
		ruleNumber
		{ after(grammarAccess.getShapeElementAccess().getXNumberParserRuleCall_0_2_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeElement__YAssignment_1_2
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getShapeElementAccess().getYNumberParserRuleCall_1_2_0()); }
		ruleNumber
		{ after(grammarAccess.getShapeElementAccess().getYNumberParserRuleCall_1_2_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeElement__WidthAssignment_2_2
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getShapeElementAccess().getWidthNumberParserRuleCall_2_2_0()); }
		ruleNumber
		{ after(grammarAccess.getShapeElementAccess().getWidthNumberParserRuleCall_2_2_0()); }
	)
;
finally {
	restoreStackSize(stackSize);
}

rule__ShapeElement__HeightAssignment_3_2
	@init {
		int stackSize = keepStackSize();
	}
:
	(
		{ before(grammarAccess.getShapeElementAccess().getHeightNumberParserRuleCall_3_2_0()); }
		ruleNumber
		{ after(grammarAccess.getShapeElementAccess().getHeightNumberParserRuleCall_3_2_0()); }
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
		{ before(grammarAccess.getPropertyAccess().getValueNumberValueParserRuleCall_2_1_0()); }
		ruleNumberValue
		{ after(grammarAccess.getPropertyAccess().getValueNumberValueParserRuleCall_2_1_0()); }
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
		{ before(grammarAccess.getPropertyAccess().getValueBooleanValueParserRuleCall_2_2_0()); }
		ruleBooleanValue
		{ after(grammarAccess.getPropertyAccess().getValueBooleanValueParserRuleCall_2_2_0()); }
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
