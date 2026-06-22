/*******************************************************************************
 * Copyright (c) 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0  
 *******************************************************************************/
package org.eclipse.elk.core.data;

import java.util.Map;
import java.util.function.Consumer;

import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.options.PortLabelPlacement;
import org.eclipse.elk.core.options.SizeOptions;
import org.eclipse.elk.core.util.IGraphElementVisitor;
import org.eclipse.elk.graph.ElkGraphElement;
import org.eclipse.elk.graph.properties.IProperty;

/**
 * Replaces (and removes) any deprecated layout options from {@link CoreOptions} with a corresponding new layout option
 * (or whatever is required to reconstruct the old behavior).
 */
public class DeprecatedLayoutOptionReplacer implements IGraphElementVisitor {

    /**
     * Rule to replace {@link CoreOptions#PORT_LABELS_NEXT_TO_PORT_IF_POSSIBLE}.
     */
    @SuppressWarnings("deprecation")
    private static final Consumer<ElkGraphElement> NEXT_TO_PORT_IF_POSSIBLE = (e) -> {
        e.getProperty(CoreOptions.PORT_LABELS_PLACEMENT).add(PortLabelPlacement.NEXT_TO_PORT_IF_POSSIBLE);
        e.setProperty(CoreOptions.PORT_LABELS_NEXT_TO_PORT_IF_POSSIBLE, null);
    };

    /**
     * Rule to move the deprecated enum value {@link SizeOptions#SPACE_EFFICIENT_PORT_LABELS} to
     * {@link PortLabelPlacement}.
     */
    @SuppressWarnings("deprecation")
    private static final Consumer<ElkGraphElement> SPACE_EFFICIENT = (e) -> {
        if (e.getProperty(CoreOptions.NODE_SIZE_OPTIONS).contains(SizeOptions.SPACE_EFFICIENT_PORT_LABELS)) {
            e.getProperty(CoreOptions.PORT_LABELS_PLACEMENT).add(PortLabelPlacement.SPACE_EFFICIENT);
            e.getProperty(CoreOptions.NODE_SIZE_OPTIONS).remove(SizeOptions.SPACE_EFFICIENT_PORT_LABELS);
        }
    };

    /**
     * Mapping of deprecated layout options to replacing rules.
     */
    @SuppressWarnings("deprecation")
    public static final Map<IProperty<?>, Consumer<ElkGraphElement>> RULES =
            Map.of(
                CoreOptions.PORT_LABELS_NEXT_TO_PORT_IF_POSSIBLE, NEXT_TO_PORT_IF_POSSIBLE,
                CoreOptions.NODE_SIZE_OPTIONS, SPACE_EFFICIENT
            );

    @Override
    public void visit(final ElkGraphElement element) {
        RULES.forEach((option, replacer) -> {
            if (element.hasProperty(option)) {
                replacer.accept(element);
                // Note that the option is _not_ removed here. The replacer has to take care of it;
                // It depends on the concrete layout option whether it should be removed or not.
            }
        });
    }

}
