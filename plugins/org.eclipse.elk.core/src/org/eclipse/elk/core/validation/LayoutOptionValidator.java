/*******************************************************************************
 * Copyright (c) 2016 TypeFox GmbH and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.validation;

import java.util.ArrayList;
import java.util.Collections;
import java.util.EnumSet;
import java.util.List;
import java.util.Map;

import org.eclipse.elk.core.data.LayoutOptionData;
import org.eclipse.elk.graph.ElkGraphElement;
import org.eclipse.elk.graph.properties.IProperty;
import org.eclipse.elk.graph.properties.IPropertyValueProxy;

/**
 * A validator for layout options that checks lower and upper bounds as well as other constraints.
 */
public class LayoutOptionValidator implements IValidatingGraphElementVisitor {
    
    /** The issues found by this validator. */
    private final List<GraphIssue> issues = new ArrayList<GraphIssue>();

    /**
     * Check the lower and upper bounds of all properties attached to graph elements.
     */
    @SuppressWarnings("unchecked")
    @Override
    public void visit(final ElkGraphElement element) {
        for (Map.Entry<IProperty<?>, Object> entry : element.getProperties()) {
            IProperty<?> property = entry.getKey();
            if (property != null) {
                Object value = entry.getValue();
                if (value instanceof IPropertyValueProxy) {
                    value = ((IPropertyValueProxy) value).resolveValue(property);
                    if (value != null) {
                        entry.setValue(value);
                    }
                }
                issues.addAll(checkProperty((IProperty<Object>) property, value, element));
            }
        }
    }
    
    /**
     * Check the lower and upper bounds of the given property.
     */
    public List<GraphIssue> checkProperty(final IProperty<Object> property, final Object value,
            final ElkGraphElement element) {
        String optionName = null;
        if (property instanceof LayoutOptionData) {
            LayoutOptionData optionData = (LayoutOptionData) property;
            optionName = optionData.getName();
            if (value != null) {
                if (!isValidType(optionData, value)) {
                    String message = "The assigned value " + value.toString()
                            + " of the option '" + optionName
                            + "' does not match the type " + optionData.getOptionClass().getSimpleName() + ".";
                    return Collections.singletonList(new GraphIssue(element, message, GraphIssue.Severity.ERROR));
                }
            }
        }
        if (value != null) {
            if (property.getLowerBound().compareTo(value) > 0) {
                if (optionName == null) {
                    optionName = property.getId();
                }
                String message = "The assigned value " + value.toString()
                        + " of the option '" + optionName
                        + "' is less than the lower bound " + property.getLowerBound().toString() + ".";
                return Collections.singletonList(new GraphIssue(element, message, GraphIssue.Severity.ERROR));
            } else if (property.getUpperBound().compareTo(value) < 0) {
                if (optionName == null) {
                    optionName = property.getId();
                }
                String message = "The assigned value " + value.toString()
                        + " of the option '" + optionName
                        + "' is greater than the upper bound " + property.getUpperBound().toString() + ".";
                return Collections.singletonList(new GraphIssue(element, message, GraphIssue.Severity.ERROR));
            }
        }
        return Collections.emptyList();
    }
    
    /**
     * @return whether the given value has the correct type
     */
    protected boolean isValidType(final LayoutOptionData optionData, final Object value) {
        switch (optionData.getType()) {
            case STRING:
                return value instanceof String;
            case BOOLEAN:
                return value instanceof Boolean;
            case INT:
                return value instanceof Integer;
            case DOUBLE:
                return value instanceof Double;
            case ENUMSET:
                return value instanceof EnumSet<?>;
            default:
                return optionData.getOptionClass().isInstance(value);
        }
    }
    
    @Override
    public List<GraphIssue> getIssues() {
        return issues;
    }

}
