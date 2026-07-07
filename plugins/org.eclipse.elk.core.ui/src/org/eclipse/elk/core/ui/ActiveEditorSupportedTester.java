/*******************************************************************************
 * Copyright (c) 2010, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.ui;

import org.eclipse.core.expressions.PropertyTester;
import org.eclipse.elk.core.service.LayoutConnectorsService;
import org.eclipse.ui.IWorkbenchPart;

/**
 * Tester class for checking whether the active editor is supported by ELK.
 * 
 * @author jjc
 */
public class ActiveEditorSupportedTester extends PropertyTester {

    @Override
    public boolean test(final Object receiver, final String property,
            final Object[] args, final Object expectedValue) {
        if (receiver instanceof IWorkbenchPart) {
            IWorkbenchPart workbenchPart = (IWorkbenchPart) receiver;
            return LayoutConnectorsService.getInstance().getConnector(workbenchPart, null) != null;
        }
        return false;
    }

}
